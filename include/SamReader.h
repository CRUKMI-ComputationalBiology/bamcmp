#ifndef SAMREADER_H
#define SAMREADER_H

#include <string>
#include <htslib/sam.h>

class SamReader
{
    public:
        SamReader(htsFile* _hf, bam_hdr_t* _header, const char* fname);
        virtual ~SamReader();
        bool is_eof() const;
        void next();
        bam1_t* getRec();
    protected:
    private:
        htsFile* hf;
        bam_hdr_t* header;
        bam1_t *prev_rec;
        bam1_t *rec;
        bool eof;
        std::string filename;
};

#endif // SAMREADER_H
