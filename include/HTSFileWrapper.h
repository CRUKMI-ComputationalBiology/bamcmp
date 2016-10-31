#ifndef HTSFILEWRAPPER_H
#define HTSFILEWRAPPER_H

#include <string>
#include <vector>
#include <htslib/hts.h>
#include <htslib/sam.h>

class HTSFileWrapper
{
    public:
        static HTSFileWrapper* begin_or_die(const char* fname, const char* mode, bam_hdr_t* header, int inputNumber, int nthreads);
        static void close(HTSFileWrapper* f);
        HTSFileWrapper(const std::string& _fname, const char* _mode, int _nthreads);
        virtual ~HTSFileWrapper();
        void checkStarted();
        void setHeader1(bam_hdr_t* h1);
        void setHeader2(bam_hdr_t* h2);
        void write1(int headerNum, bam1_t* rec);
        void ref();
        uint32_t unref();
    protected:
    private:
        static std::vector<std::pair<std::string, HTSFileWrapper*> > openOutputs;
        std::string fname;
        const char* mode;
        htsFile* hts;
        uint32_t refCount;
        int nthreads;
        int header2_offset;
        bam_hdr_t* header1;
        bam_hdr_t* header2;
        bam_hdr_t* headerOut;
        void checkHeaderNotWritten();
};

#endif // HTSFILEWRAPPER_H
