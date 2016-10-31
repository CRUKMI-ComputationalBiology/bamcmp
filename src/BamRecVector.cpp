#include "BamRecVector.h"

#include <algorithm>
#include <htslib/hts.h>

#include "util.h"

BamRecVector::BamRecVector()
{
    //ctor
}

BamRecVector::~BamRecVector()
{
    clear();
}

void BamRecVector::take_add(bam1_t* src)
{
    recs.push_back(src);
}

void BamRecVector::copy_add(bam1_t* src)
{
    recs.push_back(bam_dup1(src));
}

void BamRecVector::clear()
{
    for(std::vector<bam1_t*>::iterator it = recs.begin(), itend = recs.end(); it != itend; ++it)
    {
        bam_destroy1(*it);
    }
    recs.clear();
}

void BamRecVector::sort()
{
    std::sort(recs.begin(), recs.end(), bamrec_lt);
}

unsigned int BamRecVector::size()
{
    return recs.size();
}

bam1_t* BamRecVector::get(int index)
{
    return recs[index];
}
