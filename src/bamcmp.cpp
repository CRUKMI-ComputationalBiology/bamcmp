#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "util.h"
#include "HTSFileWrapper.h"
#include "SamReader.h"
#include "BamRecVector.h"

extern bool mixed_ordering;

enum scoringmethods
{
    scoringmethod_nmatches,
    scoringmethod_astag,
    scoringmethod_mapq,
    scoringmethod_balwayswins
};

static scoringmethods scoringmethod;

static bool warned_nm_anomaly = false;
static bool warned_nm_md_tags = false;

static void usage()
{
    fprintf(stderr, "Usage: bamcmp -1 input1.s/b/cram -2 input2.s/b/cram [-a first_only.xam] [-b second_only.xam] [-A first_better.xam] [-B second_better.xam] [-C first_worse.xam] [-D second_worse.xam] [-t nthreads] [-n | -N] [-s scoring_method]\n");
    fprintf(stderr, "\t-n\tExpect input sorted as per samtools -n (Mixed string / integer ordering, default)\n");
    fprintf(stderr, "\t-N\tExpect input sorted as per Picard / htsjdk name ordering (lexical ordering)\n");
    fprintf(stderr, "\t-s match\tScore hits by the number of bases that match the reference, as given by the CIGAR string and NM / MD attributes (default)\n");
    fprintf(stderr, "\t-s as\tScore hits according to the AS attribute written by some aligners\n");
    fprintf(stderr, "\t-s mapq\tScore hits according to the MAPQ SAM field\n");
    fprintf(stderr, "\t-s balwayswins\tAlways award hits to input B, regardless of alignment scores (equivalent to filtering A by any read mapped in B)\n");
    exit(1);
}

static bool aux_is_int(uint8_t* rec)
{
    switch(*rec)
    {
    case 'c':
    case 'C':
    case 's':
    case 'S':
    case 'i':
    case 'I':
        return true;
    default:
        return false;
    }
}

static uint32_t get_alignment_score(bam1_t* rec, bool is_input_a)
{
    switch(scoringmethod)
    {
    case scoringmethod_nmatches:
    {
        bool seen_equal_or_diff = false;
        int32_t cigar_total = 0;
        const uint32_t* cigar = bam_get_cigar(rec);
        int32_t indel_edit_distance = 0;

        for(int i = 0; i < rec->core.n_cigar; ++i)
        {
            // CIGAR scoring: score points for matching bases, and negatives for deletions
            // since otherwise 10M10D10M would score the same as 20M. Insertions, clipping etc
            // don't need to score a penalty since they skip bases in the query.
            // CREF_SKIP (N / intron-skip operator) is acceptable: 10M1000N10M is as good as 20M.
            // Insertions are counted to correct the NM tag below only.

            int32_t n = bam_cigar_oplen(cigar[i]);
            switch(bam_cigar_op(cigar[i]))
            {
            case BAM_CEQUAL:
                seen_equal_or_diff = true;
            // fall through
            case BAM_CMATCH:
                cigar_total += n;
                break;

            case BAM_CDEL:
                indel_edit_distance += n;
                cigar_total -= n;
                break;

            case BAM_CDIFF:
                seen_equal_or_diff = true;
                break;

            case BAM_CINS:
                indel_edit_distance += n;
                break;

            default:
                break;
            }
        }

        // The BAM_CMATCH operator (unlike BAM_CEQUAL or BAM_CDIFF) could mean a match or a mismatch
        // with same length (e.g. a SNP). If the file doesn't seem to use the advanced operators try to
        // spot mismatches from metadata tags.
        if(!seen_equal_or_diff)
        {
            uint8_t* nm_rec = bam_aux_get(rec, "NM");
            if(nm_rec && aux_is_int(nm_rec))
            {
                int32_t nm = bam_aux2i(nm_rec);
                if(nm < indel_edit_distance)
                {
                    if(!warned_nm_anomaly)
                    {
                        fprintf(stderr, "Warning: anomaly in record %s: NM is %d but there are at least %d indel bases in the CIGAR string\n", bam_get_qname(rec), nm, indel_edit_distance);
                        fprintf(stderr, "There may be more records with this problem, but the warning will not be repeated\n");
                        warned_nm_anomaly = true;
                    }
                }
                else
                {
                    cigar_total -= (bam_aux2i(nm_rec) - indel_edit_distance);
                    seen_equal_or_diff = true;
                }
            }
        }

        if(!seen_equal_or_diff)
        {
            uint8_t* md_rec = bam_aux_get(rec, "MD");
            if(md_rec)
            {
                char* mdstr = bam_aux2Z(md_rec);
                if(mdstr)
                {
                    seen_equal_or_diff = true;
                    bool in_deletion = false;
                    for(; *mdstr; ++mdstr)
                    {
                        // Skip deletions, which are already penalised.
                        // Syntax seems to be: numbers mean base strings that match the reference; ^ followed by letters means
                        // a deletion; letters without the preceding ^ indicate a mismatch.
                        char c = *mdstr;
                        if(c == '^')
                        {
                            in_deletion = true;
                        }
                        else if(isdigit(c))
                        {
                            in_deletion = false;
                        }
                        else if(!in_deletion)
                        {
                            // Mismatch
                            cigar_total--;
                        }
                    }
                }
            }
        }

        if((!seen_equal_or_diff) && !warned_nm_md_tags)
        {
            fprintf(stderr, "Warning: input file does not use the =/X CIGAR operators, or include NM or MD tags, so I have no way to spot length-preserving reference mismatches.\n");
            fprintf(stderr, "At least record %s exhibited this problem; there may be others but the warning will not be repeated. I will assume M CIGAR operators indicate a match.\n", bam_get_qname(rec));
            warned_nm_md_tags = true;
        }
        return std::max(cigar_total, 0);
    }
    // End the CIGAR string scoring method. Thankfully the others are much simpler to implement:

    case scoringmethod_astag:
    {
        uint8_t* score_rec = bam_aux_get(rec, "AS");
        if(!score_rec)
        {
            fprintf(stderr, "Fatal: At least record %s doesn't have an AS tag as required.\n", bam_get_qname(rec));
            exit(1);
        }
        return bam_aux2i(score_rec);
    }

    case scoringmethod_mapq:
        return rec->core.qual;

    case scoringmethod_balwayswins:
        // Mapped B records beat any A record, beats an unmapped B record.
        if(is_input_a)
        {
            return 1;
        }
        else if(!(rec->core.flag & BAM_FUNMAP))
        {
            return 2;
        }
        else
        {
            return 0;
        }
    }
    return -1;
}

static bool uniqueValue(const std::vector<HTSFileWrapper*>& in)
{

    bool outValid = false;
    HTSFileWrapper* out = 0;

    for(std::vector<HTSFileWrapper*>::const_iterator it = in.begin(), itend = in.end(); it != itend; ++it)
    {
        if(!outValid)
        {
            out = *it;
            outValid = true;
        }
        else if(out != *it)
        {
            return false;
        }
    }
    return outValid;
}

static void clearMateInfo(BamRecVector& v)
{

    for(int i = 0, ilim = v.size(); i != ilim; ++i)
    {
        uint32_t maten = (uint32_t)flag2mate(v.get(i));
        bam_aux_append(v.get(i), "om", 'i', sizeof(uint32_t), (uint8_t*)&maten);

        v.get(i)->core.flag &= ~(BAM_FPROPER_PAIR | BAM_FMREVERSE | BAM_FPAIRED | BAM_FMUNMAP | BAM_FREAD1 | BAM_FREAD2);
        v.get(i)->core.mtid = -1;
        v.get(i)->core.mpos = -1;
    }
}

int main(int argc, char** argv)
{
    char *in1_name = NULL, *in2_name = NULL, *firstbetter_name = NULL, *secondbetter_name = NULL,
          *firstworse_name = NULL, *secondworse_name = NULL, *first_name = NULL, *second_name = NULL;

    int nthreads = 1;
    std::string scoring_method_string = "match";

    int c;
    while ((c = getopt(argc, argv, "a:b:1:2:t:A:B:C:D:nNs:")) >= 0)
    {
        switch (c)
        {
        case '1':
            in1_name = optarg;
            break;
        case '2':
            in2_name = optarg;
            break;
        case 'a':
            first_name = optarg;
            break;
        case 'b':
            second_name = optarg;
            break;
        case 'A':
            firstbetter_name = optarg;
            break;
        case 'B':
            secondbetter_name = optarg;
            break;
        case 'C':
            firstworse_name = optarg;
            break;
        case 'D':
            secondworse_name = optarg;
            break;
        case 't':
            nthreads = atoi(optarg);
            break;
        case 'n':
            mixed_ordering = true;
            break;
        case 'N':
            mixed_ordering = false;
            break;
        case 's':
            scoring_method_string = std::string(optarg);
            break;
        default:
            usage();
        }
    }

    if(in1_name == NULL)
    {
        usage();
    }
    if(in2_name == NULL)
    {
        usage();
    }
    if(!(first_name != NULL || second_name != NULL || firstbetter_name != NULL || secondbetter_name != NULL))
    {
        fprintf(stderr, "bamcmp is useless without at least one of -1, -2, -A or -B\n");
        usage();
    }

    if(scoring_method_string.compare("match") == 0)
    {
        scoringmethod = scoringmethod_nmatches;
    }
    else if(scoring_method_string.compare("mapq") == 0)
    {
        scoringmethod = scoringmethod_mapq;
    }
    else if(scoring_method_string.compare("as") == 0)
    {
        scoringmethod = scoringmethod_astag;
    }
    else if(scoring_method_string.compare("balwayswins") == 0)
    {
        scoringmethod = scoringmethod_balwayswins;
    }
    else
    {
        usage();
    }

    htsFile *in1hf = NULL, *in2hf = NULL;
    HTSFileWrapper *firstbetter_out = NULL, *secondbetter_out = NULL, *firstworse_out = NULL, *secondworse_out = NULL, *first_out = NULL, *second_out = NULL;
    in1hf = hts_begin_or_die(in1_name, "r", 0, nthreads);
    in2hf = hts_begin_or_die(in2_name, "r", 0, nthreads);

    bam_hdr_t* header1 = sam_hdr_read(in1hf);
    bam_hdr_t* header2 = sam_hdr_read(in2hf);

    // Permit the outputs using like headers to share a file if they gave the same name.

    if(firstbetter_name)
    {
        firstbetter_out = HTSFileWrapper::begin_or_die(firstbetter_name, "wb0", header1, 1, nthreads);
    }
    if(secondbetter_name)
    {
        secondbetter_out = HTSFileWrapper::begin_or_die(secondbetter_name, "wb0", header2, 2, nthreads);
    }
    if(firstworse_name)
    {
        firstworse_out = HTSFileWrapper::begin_or_die(firstworse_name, "wb0", header1, 1, nthreads);
    }
    if(secondworse_name)
    {
        secondworse_out = HTSFileWrapper::begin_or_die(secondworse_name, "wb0", header2, 2, nthreads);
    }
    if(first_name)
    {
        first_out = HTSFileWrapper::begin_or_die(first_name, "wb0", header1, 1, nthreads);
    }
    if(second_name)
    {
        second_out = HTSFileWrapper::begin_or_die(second_name, "wb0", header2, 2, nthreads);
    }

    SamReader in1(in1hf, header1, in1_name);
    SamReader in2(in2hf, header2, in2_name);

    BamRecVector seqs1, seqs2;
    std::vector<HTSFileWrapper*> seqs1Files, seqs2Files;

    while((!in1.is_eof()) && (!in2.is_eof()))
    {
        std::string qname1 = std::string(bam_get_qname(in1.getRec()));
        std::string qname2 = std::string(bam_get_qname(in2.getRec()));

        if(qname1 == qname2)
        {
            seqs1.clear();
            seqs2.clear();
            seqs1Files.clear();
            seqs2Files.clear();

            std::string qn;
            while((!in1.is_eof()) && (qn = bam_get_qname(in1.getRec())) == qname1)
            {
                seqs1.copy_add(in1.getRec());
                in1.next();
            }

            seqs1.sort();
            seqs1Files.resize(seqs1.size(), 0);

            while((!in2.is_eof()) && (qn = bam_get_qname(in2.getRec())) == qname2)
            {
                seqs2.copy_add(in2.getRec());
                in2.next();
            }

            seqs2.sort();
            seqs2Files.resize(seqs2.size(), 0);

            unsigned int idx1 = 0, idx2 = 0;
            while(idx1 < seqs1.size() && idx2 < seqs2.size())
            {
                if(bamrec_eq(seqs1.get(idx1), seqs2.get(idx2)))
                {
                    uint32_t score1 = 0;
                    uint32_t score2 = 0;

                    int group_start_idx1 = idx1, group_start_idx2 = idx2;

                    score1 = get_alignment_score(seqs1.get(idx1), true);
                    score2 = get_alignment_score(seqs2.get(idx2), false);

                    // Either input may have multiple candidate matches. Compare the best match found in each group
                    // and then emit the whole group as firstbetter or secondbetter.

                    while(idx1 + 1 < seqs1.size() && bamrec_eq(seqs1.get(group_start_idx1), seqs1.get(idx1 + 1)))
                    {
                        ++idx1;
                        score1 = std::max(score1, get_alignment_score(seqs1.get(idx1), true));
                    }

                    while(idx2 + 1 < seqs2.size() && bamrec_eq(seqs1.get(group_start_idx1), seqs2.get(idx2 + 1)))
                    {
                        ++idx2;
                        score2 = std::max(score2, get_alignment_score(seqs2.get(idx2), false));
                    }

                    for(uint32_t i = group_start_idx1; i <= idx1; ++i)
                    {
                        bam_aux_append(seqs1.get(i), "as", 'i', sizeof(uint32_t), (uint8_t*)&score1);
                        bam_aux_append(seqs1.get(i), "bs", 'i', sizeof(uint32_t), (uint8_t*)&score2);
                    }

                    for(uint32_t i = group_start_idx2; i <= idx2; ++i)
                    {
                        bam_aux_append(seqs2.get(i), "as", 'i', sizeof(uint32_t), (uint8_t*)&score1);
                        bam_aux_append(seqs2.get(i), "bs", 'i', sizeof(uint32_t), (uint8_t*)&score2);
                    }

                    HTSFileWrapper *firstRecordsFile, *secondRecordsFile;

                    if(score1 > score2)
                    {
                        firstRecordsFile = firstbetter_out;
                        secondRecordsFile = secondworse_out;
                    }
                    else
                    {
                        firstRecordsFile = firstworse_out;
                        secondRecordsFile = secondbetter_out;
                    }

                    for(uint32_t i = group_start_idx1; i <= idx1; ++i)
                    {
                        seqs1Files[i] = firstRecordsFile;
                    }

                    for(uint32_t i = group_start_idx2; i <= idx2; ++i)
                    {
                        seqs2Files[i] = secondRecordsFile;
                    }
                    ++idx1;
                    ++idx2;
                }
                else if(bamrec_lt(seqs1.get(idx1), seqs2.get(idx2)))
                {
                    seqs1Files[idx1] = first_out;
                    ++idx1;
                }
                else
                {
                    seqs2Files[idx2] = second_out;
                    ++idx2;
                }
            }

            for(; idx1 < seqs1.size(); ++idx1)
            {
                seqs1Files[idx1] = first_out;
            }

            for(; idx2 < seqs2.size(); ++idx2)
            {
                seqs2Files[idx2] = second_out;
            }

            // Figure out whether we're splitting the mates up in either case.
            // If they are split up, clear mate information to make the file consistent.
            if(!uniqueValue(seqs1Files))
            {
                clearMateInfo(seqs1);
            }

            if(!uniqueValue(seqs2Files))
            {
                clearMateInfo(seqs2);
            }

            for(int i = 0, ilim = seqs1.size(); i != ilim; ++i)
            {
                if(seqs1Files[i])
                {
                    seqs1Files[i]->write1(1, seqs1.get(i));
                }
            }
            for(int i = 0, ilim = seqs2.size(); i != ilim; ++i)
            {
                if(seqs2Files[i])
                {
                    seqs2Files[i]->write1(2, seqs2.get(i));
                }
            }
        }
        else if(qname_cmp(qname1.c_str(), qname2.c_str()) < 0)
        {
            if(first_out)
            {
                first_out->write1(1, in1.getRec());
            }
            in1.next();
        }
        else
        {
            if(second_out)
            {
                second_out->write1(2, in2.getRec());
            }
            in2.next();
        }
    }

    // One or other file has reached EOF. Write the remainder as first- or second-only records.
    if(first_out)
    {
        while(!in1.is_eof())
        {
            first_out->write1(1, in1.getRec());
            in1.next();
        }
    }

    if(second_out)
    {
        while(!in2.is_eof())
        {
            second_out->write1(2, in2.getRec());
            in2.next();
        }
    }

    hts_close(in1hf);
    hts_close(in2hf);
    if(first_out)
    {
        HTSFileWrapper::close(first_out);
    }
    if(second_out)
    {
        HTSFileWrapper::close(second_out);
    }
    if(firstbetter_out)
    {
        HTSFileWrapper::close(firstbetter_out);
    }
    if(secondbetter_out)
    {
        HTSFileWrapper::close(secondbetter_out);
    }
    if(firstworse_out)
    {
        HTSFileWrapper::close(firstworse_out);
    }
    if(secondworse_out)
    {
        HTSFileWrapper::close(secondworse_out);
    }
}
