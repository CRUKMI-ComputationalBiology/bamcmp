#include "HTSFileWrapper.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

std::vector<std::pair<std::string, HTSFileWrapper*> > HTSFileWrapper::openOutputs;

HTSFileWrapper* HTSFileWrapper::begin_or_die(const char* fname, const char* mode, bam_hdr_t* header, int inputNumber, int nthreads)
{
    if(inputNumber != 1 && inputNumber != 2)
    {
        fprintf(stderr, "inputNumber must be 1 or 2\n");
        exit(1);
    }
    HTSFileWrapper* ret = NULL;
    std::string sfname(fname);
    for(std::vector<std::pair<std::string, HTSFileWrapper*> >::iterator it = HTSFileWrapper::openOutputs.begin(), itend = HTSFileWrapper::openOutputs.end(); it != itend && !ret; ++it)
    {
        if(it->first == sfname)
        {
            it->second->ref();
            ret = it->second;
        }
    }
    if(!ret)
    {
        ret = new HTSFileWrapper(sfname, mode, nthreads);
        HTSFileWrapper::openOutputs.push_back(std::make_pair(sfname, ret));
    }
    if(inputNumber == 1)
    {
        ret->setHeader1(header);
    }
    else
    {
        ret->setHeader2(header);
    }
    return ret;
}

void HTSFileWrapper::close(HTSFileWrapper* f)
{
    if(!f->unref())
    {
        delete f;
    }
}

HTSFileWrapper::HTSFileWrapper(const std::string& _fname, const char* _mode, int _nthreads)  :
        fname(_fname), mode(_mode), hts(0), refCount(1), nthreads(_nthreads), header2_offset(0), header1(0), header2(0), headerOut(0)
{
    //ctor
}

HTSFileWrapper::~HTSFileWrapper()
{
    //dtor
}

void HTSFileWrapper::setHeader1(bam_hdr_t* h1)
{
    checkHeaderNotWritten();
    header1 = h1;
}

void HTSFileWrapper::setHeader2(bam_hdr_t* h2)
{
    checkHeaderNotWritten();
    header2 = h2;
}

void HTSFileWrapper::ref()
{
    ++refCount;
}

uint32_t HTSFileWrapper:: unref()
{
    checkStarted();
    if(refCount == 1)
    {
        hts_close(hts);
    }
    return --refCount;
}

void HTSFileWrapper::checkStarted()
{
    if(hts)
    {
        return;
    }
    if(!(header1 || header2))
    {
        fprintf(stderr, "Started writing records without any header\n");
        exit(1);
    }
    else if(!header2)
    {
        headerOut = header1;
    }
    else if(!header1)
    {
        headerOut = header2;
    }
    else
    {
        header2_offset = header1->n_targets;
        headerOut = bam_hdr_dup(header1);
        headerOut->n_targets += header2->n_targets;
        headerOut->target_len = (uint32_t*)realloc(headerOut->target_len, sizeof(uint32_t) * headerOut->n_targets);
        headerOut->target_name = (char**)realloc(headerOut->target_name, sizeof(char*) * headerOut->n_targets);

        if((!headerOut->target_len) || (!headerOut->target_name))
        {
            goto oom;
        }

        // Prefix header1 names with A_
        for(int i = 0; i < header1->n_targets; ++i)
        {
            int len = strlen(headerOut->target_name[i]);
            headerOut->target_name[i] = (char*)realloc(headerOut->target_name[i], len + 3); // A_ prefix plus null terminator
            if(!headerOut->target_name[i])
            {
                goto oom;
            }
            memmove(headerOut->target_name[i] + 2, headerOut->target_name[i], len + 1);
            memcpy(headerOut->target_name[i], "A_", 2);
        }

        // Copy header2 names; add B_ prefix
        for(int i = 0; i < header2->n_targets; ++i)
        {
            headerOut->target_len[i + header2_offset] = header2->target_len[i];
            headerOut->target_name[i + header2_offset] = (char*)malloc(strlen(header2->target_name[i]) + 3);
            if(!headerOut->target_name[i + header2_offset])
            {
                goto oom;
            }
            sprintf(headerOut->target_name[i + header2_offset], "B_%s", header2->target_name[i]);
        }

        // Find all text @SQ records in header2:
        std::vector<std::pair<int, int> > sqheaders;

        int offset = 0;
        char* sq_start;
        while((sq_start = strstr(header2->text + offset, "@SQ")) != 0)
        {
            // Must occur at the start of the header or after a newline:
            if(sq_start != header2->text && (*(sq_start - 1)) != '\n')
            {
                ++offset;
                continue;
            }

            int start_off = sq_start - header2->text;
            int end_off;

            char* nl = strchr(sq_start, '\n');
            if(!nl)
            {
                end_off = strlen(header2->text);
            }
            else
            {
                end_off = nl - header2->text;
            }

            sqheaders.push_back(std::make_pair(start_off, end_off));
            offset = end_off;
        }

        // Count how much extra we're going to add to header1:
        int additional_bytes = (2 * header1->n_targets); // A_ prefix for all existing @SQs
        for(int i = 0, ilim = sqheaders.size(); i != ilim; ++i)
        {
            additional_bytes += ((sqheaders[i].second - sqheaders[i].first) + 3); // Existing header, plus a newline, plus B_ prefix.
        }
        free(headerOut->text);
        headerOut->text = (char*)malloc(headerOut->l_text + additional_bytes + 1); // Null terminator as well?
        headerOut->l_text += additional_bytes;
        if(!headerOut->text)
        {
            goto oom;
        }

        // Rewrite: copy header1 lines, except for @SQs which gain an A_ prefix, and insert the (prefixed) header2 @SQs
        // after the last one from header1.
        char* nl;
        char* read_offset = header1->text;
        char* write_offset = headerOut->text;
        bool more_to_read = true;
        int sqs_remaining = header1->n_targets;
        while(more_to_read)
        {
            nl = strchr(read_offset, '\n');
            if(!nl)
            {
                nl = header1->text + header1->l_text;
                more_to_read = false;
            }
            else
            {
                // Include the newline when copying.
                ++nl;
            }
            if(strncmp(read_offset, "@SQ", 3) == 0)
            {
                char* sn_offset = strstr(read_offset, "SN:");
                if(sn_offset)
                {
                    // Prefix A_ to the @SQ line
                    char* insert_offset = sn_offset + 3;
                    int len = insert_offset - read_offset;
                    memcpy(write_offset, read_offset, len);

                    write_offset += len;
                    read_offset += len;

                    memcpy(write_offset, "A_", 2);
                    write_offset += 2;

                    len = nl - read_offset;
                    memcpy(write_offset, read_offset, len);

                    write_offset += len;
                    read_offset += len;

                    if(!more_to_read)
                    {
                        // Adding at the end; add newline.
                        *(write_offset++) = '\n';
                    }
                    if(!(--sqs_remaining))
                    {
                        // That was the last of header1's @SQ lines. Insert header2's lines here.
                        for(int i = 0, ilim = sqheaders.size(); i != ilim; ++i)
                        {
                            int h2readoff = sqheaders[i].first;
                            int h2readlim = sqheaders[i].second;
                            char* h2sn = strstr(header2->text + h2readoff, "SN:");
                            if(h2sn)
                            {
                                // Skip over the SN: tag
                                h2sn += 3;
                                int h2insoff = h2sn - header2->text;
                                int len = h2insoff - h2readoff;

                                memcpy(write_offset, header2->text + h2readoff, len);
                                write_offset += len;

                                memcpy(write_offset, "B_", 2);
                                write_offset += 2;

                                len = h2readlim - h2insoff;
                                memcpy(write_offset, header2->text + h2insoff, len);
                                write_offset += len;
                            }
                            else
                            {
                                int len = h2readlim - h2readoff;
                                memcpy(write_offset, header2->text + h2readoff, len);
                                write_offset += len;
                            }
                            if(more_to_read || i != ilim - 1)
                            {
                                *(write_offset++) = '\n';
                            }
                        }
                    }
                    continue;
                }
                // Else fall through to just copy the whole line:
            }
            // Just a regular header line: copy it all:
            int len = nl - read_offset;
            memcpy(write_offset, read_offset, len);
            write_offset += len;
            read_offset += len;
        }
        headerOut->text[headerOut->l_text] = '\0';
    }
    // Header complete, now open and write it:
    hts = hts_begin_or_die(fname.c_str(), mode, headerOut, nthreads);
    return;
oom:
    fprintf(stderr, "Malloc failure while building combined header\n");
    exit(1);
}

void HTSFileWrapper::write1(int headerNum, bam1_t* rec)
{
    checkStarted();
    if(headerNum == 2)
    {
        if(rec->core.tid != -1)
        {
            rec->core.tid += header2_offset;
        }
        if(rec->core.mtid != -1)
        {
            rec->core.mtid += header2_offset;
        }
    }
    if (sam_write1(hts, headerOut, rec) == -1)
    {
        // FIXME - Handle this error condition!
    }
    if(headerNum == 2)
    {
        if(rec->core.tid != -1)
        {
            rec->core.tid -= header2_offset;
        }
        if(rec->core.mtid != -1)
        {
            rec->core.mtid -= header2_offset;
        }
    }
}

void HTSFileWrapper::checkHeaderNotWritten()
{
    if(headerOut)
    {
        fprintf(stderr, "htsFileWrapper header changed after already being written\n");
        exit(1);
    }
}
