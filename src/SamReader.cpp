#include "SamReader.h"

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>

#include "util.h"

extern bool mixed_ordering;

SamReader::SamReader(htsFile* _hf, bam_hdr_t* _header, const char* fname) : hf(_hf), header(_header), eof(false), filename(fname)
{
    rec = bam_init1();
    prev_rec = bam_init1();
    next();
}

SamReader::~SamReader()
{
    bam_destroy1(rec);
}

bool SamReader::is_eof() const
{
    return eof;
}

void SamReader::next()
{
    if(eof)
    {
        return;
    }

    if(rec->data)
    {
        bam_copy1(prev_rec, rec);
    }

    if(sam_read1(hf, header, rec) < 0)
    {
        eof = true;
    }
    if(prev_rec->data && (!eof) && qname_cmp(bam_get_qname(rec), bam_get_qname(prev_rec)) < 0)
    {
        fprintf(stderr, "Order went backwards! In file %s, record %s belongs before %s. Re-sort your files and try again.\n", filename.c_str(), bam_get_qname(rec), bam_get_qname(prev_rec));
        if(mixed_ordering)
        {
            fprintf(stderr, "Expected order was the mixed string/integer ordering produced by samtools sort -n; use -N to switch to Picard / htsjdk string ordering\n");
        }
        else
        {
            fprintf(stderr, "Expected order was Picard / htsjdk string ordering; use -n to switch to samtools sort -n ordering\n");
        }
        exit(1);
    }
}

bam1_t* SamReader::getRec()
{
    return rec;
}
