#ifndef BROWHEEL_H
#define BROWHEEL_H

#include <bits/stdc++.h>
#include <limits>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <tuple>
#include <map>
#include <cctype>
#include <algorithm>
#include <filesystem>

#define M64 1
#include "../external/gsacak.h"
#include "../external/gsacak.c"
#undef max

const static std::string int2base="#$NACGT";
static std::map<char, int8_t> base2int={{'#',0}, {'$',1}, {'N',2}, {'A',3}, {'C',4}, {'G',5}, {'T',6}, {'n',2}, {'a',3}, {'c',4}, {'g',5}, {'t',6}};

void gsufsortSA(std::string SAfile, uint8_t *revref, uint64_t revref_sz)
{
    uint64_t *SA = (uint64_t *) malloc(revref_sz*sizeof(uint64_t));
    sacak(revref, SA, revref_sz);
    std::ofstream fSA(SAfile);
    fSA.write((char *)SA, revref_sz*sizeof(uint64_t));
    delete[] SA;
}

uint8_t* SA2bwt(std::string SAfile, uint64_t batch, uint8_t *revseq)
{
    uint64_t revseq_sz = std::filesystem::file_size(SAfile)/sizeof(uint64_t);
    std::ifstream fin(SAfile);
    uint64_t *SA = (uint64_t *) malloc(batch*sizeof(uint64_t));
    uint8_t *bwt = (uint8_t *) malloc(revseq_sz);
    uint64_t bwt_sz = 0;
    while (true)
    {
        fin.read((char *)SA, batch*sizeof(uint64_t));
        uint64_t gc = fin.gcount()/sizeof(uint64_t);
        if (gc==0)
            break;
        for (uint64_t i=0; i<gc; ++i, ++bwt_sz)
        {
            uint64_t j = (SA[i]>0 ? (SA[i]-1) : (revseq_sz-1));
            bwt[bwt_sz] = revseq[j];
        }
    }
    return bwt;
}

struct RankVec
{
    uint64_t** R1=NULL;
    uint16_t** R2=NULL;
    const static uint64_t b1=8;
    const static uint64_t b2=256*32;
    const static uint8_t sigma=6;
    std::array<uint64_t,sigma> C;
    uint8_t * bwt=NULL;
    uint64_t bwt_sz;

    RankVec()
    {
    }

    ~RankVec()
    {
        if (bwt)
            delete[] bwt;
        release(R1);
        release(R2);
    }

    template <typename U>
    void release(U** R)
    {
        if(R)
        {
            delete[] R[0];
            delete[] R;
        }
    }

    template <typename U>
    void apply(U** & R, uint64_t block_size)
    {
        uint64_t sz=(bwt_sz-1)/block_size+2;
        release(R);
        R=new U*[sigma];
        R[0]=new U[sigma*sz];
        for(uint8_t c=1; c<sigma; ++c)
            R[c]=R[c-1]+sz;
    }

    void count()
    {
        apply(R1, b1*b2);
        apply(R2, b1);
        for(uint8_t c=0; c<sigma; ++c)
        {
            R1[c][0]=0;
            R2[c][0]=0;
        }
        uint64_t j1=0;
        for(uint64_t i=0, j2=0; i<bwt_sz; ++i)
        {
            if(i%(b1*b2)==0)
            {
                ++j1;
                for(uint8_t c=0; c<sigma; ++c)
                    R1[c][j1]=R1[c][j1-1];
            }
            if(i%b1==0)
            {
                ++j2;
                for(uint8_t c=0; c<sigma; ++c)
                    if(j2%b2==0)
                        R2[c][j2]=0;
                    else
                        R2[c][j2]=R2[c][j2-1];
            }
            uint8_t c=bwt[i];
            if(c>0)
            {
                ++R1[c-1][j1];
                if(j2%b2!=0)
                    ++R2[c-1][j2];
            }
        }

        C[0]=0;
        for(int c=1; c<sigma; ++c)
            C[c]=C[c-1]+R1[c-1][j1];
    }

    uint64_t rank(uint8_t c, int64_t i)
    {
        if (i<0)
            return 0;
        uint64_t r=R1[c-1][i/(b1*b2)]+R2[c-1][i/b1];
        for(uint64_t j=(i/b1)*b1; j<=i; ++j)
            if(bwt[j]==c)
                ++r;
        return r;
    }

    void PreRange(int64_t& i, int64_t& j, uint8_t c)
    {
        i=C[c-1]+rank(c,i-1)+1;
        j=C[c-1]+rank(c,j);
    }

    void saveRankVec(std::string bwtFile, std::string bwtRankFile)
    {
        std::ofstream fout(bwtFile);
        fout.write((char*)bwt, bwt_sz);
        fout.close();
        
        fout.open(bwtRankFile);
        fout.write((char*)C.data(), C.size() * sizeof(C[0]));
        uint64_t sz=(bwt_sz-1)/(b1*b2)+2;
        fout.write((char*)R1[0],this->sigma*sz*sizeof(R1[0][0]));
        sz=(bwt_sz-1)/b1+2;
        fout.write((char*)R2[0],this->sigma*sz*sizeof(R2[0][0]));
    }

    void loadRankVec(std::string bwtFile, std::string bwtRankFile)
    {
        bwt_sz = std::filesystem::file_size(bwtFile);
        bwt = (uint8_t *) malloc(bwt_sz);
        std::ifstream fin(bwtFile);
        fin.read((char*)bwt, bwt_sz);
        fin.close();

        apply(R1, b1*b2);
        apply(R2, b1);
        fin.open(bwtRankFile);
        fin.read((char*)C.data(), C.size() * sizeof(C[0]));
        uint64_t sz=(bwt_sz+b1*b2-1)/(b1*b2)+1;
        fin.read((char*)R1[0],this->sigma*sz*sizeof(R1[0][0]));
        sz=(bwt_sz-1)/b1+2;
        fin.read((char*)R2[0],this->sigma*sz*sizeof(R2[0][0]));
    }
};


bool check_bwt(RankVec &bwtRank, uint8_t *revref)
{
    for(int64_t i=0, j=bwtRank.bwt_sz-2; j>=0; --j)
    {
        if (bwtRank.bwt[i]!=revref[j])
            return false;
        int c=bwtRank.bwt[i];
        i=bwtRank.C[c-1]+bwtRank.rank(c,i);
    }
    return true;
}

// int64_t SimSuffix(int64_t i)
// {
//     int j=0;
//     while(sra[i]==0)
//     {
//         int c=bwt[i];
//         i=C[c-1]+bwtRank.rank(c,i,bwt);
//         ++j;
//     }
//     return ssa[sraRank.rank(1,i,sra)-1]+j;
// }

#endif
