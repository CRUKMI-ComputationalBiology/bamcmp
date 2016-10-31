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
