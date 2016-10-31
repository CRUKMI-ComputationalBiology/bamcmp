#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <htslib/sam.h>

htsFile* hts_begin_or_die(const char* filename, const char* mode, bam_hdr_t* header, int nthreads);
int strnum_cmp(const char *_a, const char *_b);
int qname_cmp(const char* qa, const char* qb);
int flag2mate(const bam1_t* rec);
bool bamrec_eq(const bam1_t* a, const bam1_t* b);
bool bamrec_lt(const bam1_t* a, const bam1_t* b);

#endif // UTIL_H_INCLUDED
