/*
bamcmp, A program separates human and non-human DNA or RNA from contaminated
short reads produced by NGS

Copyright (C) 2016  Christopher Smowton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


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
