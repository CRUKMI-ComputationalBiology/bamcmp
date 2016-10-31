#ifndef BAMRECVECTOR_H
#define BAMRECVECTOR_H

#include <vector>
#include <htslib/sam.h>

class BamRecVector
{
    public:
        BamRecVector();
        virtual ~BamRecVector();
        void take_add(bam1_t* src);
        void copy_add(bam1_t* src);
        void clear();
        void sort();
        unsigned int size();
        bam1_t* get (int index);
    protected:
    private:
        std::vector<bam1_t*> recs;
};

#endif // BAMRECVECTOR_H
