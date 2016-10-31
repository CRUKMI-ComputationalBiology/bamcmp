#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

bool mixed_ordering = true;

htsFile* hts_begin_or_die(const char* filename, const char* mode, bam_hdr_t* header, int nthreads)
{
    htsFile* hf = hts_open(filename, mode);
    if(hf == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", filename);
        exit(1);
    }
    // header should be 0 if we're opening to read
    if(header != NULL && sam_hdr_write(hf, header))
    {
        fprintf(stderr, "Failed to write header for %s\n", filename);
        exit(1);
    }
    if(nthreads != 1)
    {
        hts_set_threads(hf, nthreads);
    }

    return hf;
}

// Borrowed from Samtools source, since samtools sort -n uses this ordering:
int strnum_cmp(const char *_a, const char *_b)
{
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb)
    {
        if (isdigit(*pa) && isdigit(*pb))
        {
            while (*pa == '0') {++pa;}
            while (*pb == '0') {++pb;}
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb)
            {
                ++pa;
                ++pb;
            }
            if (isdigit(*pa) && isdigit(*pb))
            {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i]))
                {
                    ++i;
                }
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            }
            else if (isdigit(*pa))
            {
                return 1;
            }
            else if (isdigit(*pb))
            {
                return -1;
            }
            else if (pa - a != pb - b)
            {
                return pa - a < pb - b? 1 : -1;
            }
        }
        else
        {
            if (*pa != *pb) {return (int)*pa - (int)*pb;}
            ++pa;
            ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

int qname_cmp(const char* qa, const char* qb)
{
    if(mixed_ordering)
    {
        return strnum_cmp(qa, qb);
    }
    else
    {
        return strcmp(qa, qb);
    }
}

int flag2mate(const bam1_t* rec)
{
    if(rec->core.flag & BAM_FREAD1)
    {
        return 1;
    }
    else if(rec->core.flag & BAM_FREAD2)
    {
        return 2;
    }
    return 0;
}

bool bamrec_eq(const bam1_t* a, const bam1_t* b)
{
    return flag2mate(a) == flag2mate(b);
}

bool bamrec_lt(const bam1_t* a, const bam1_t* b)
{
    return flag2mate(a) < flag2mate(b);
}
