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

using NUCTYPE = uint8_t;
using SIZETYPE = uint64_t;
using R2TYPE = uint8_t;

const static std::string NUCTYPE2base="#$NACGT";
static std::map<char, NUCTYPE> base2NUCTYPE={{'#',0}, {'$',1}, {'N',2}, {'A',3}, {'C',4}, {'G',5}, {'T',6}, {'n',2}, {'a',3}, {'c',4}, {'g',5}, {'t',6}};

void gsufsortSA(std::string SAfile, NUCTYPE *revref, SIZETYPE revref_sz)
{
    SIZETYPE *SA = new SIZETYPE[revref_sz];
    sacak(revref, SA, revref_sz);
    std::ofstream fSA(SAfile);
    fSA.write((char *)SA, revref_sz*sizeof(SIZETYPE));
    delete[] SA;
}

NUCTYPE* SA2bwt(std::string SAfile, SIZETYPE batch, NUCTYPE *revseq)
{
    SIZETYPE revseq_sz = std::filesystem::file_size(SAfile)/sizeof(SIZETYPE);
    std::ifstream fin(SAfile);
    SIZETYPE *SA = new SIZETYPE[batch];
    NUCTYPE *bwt = new NUCTYPE[revseq_sz];
    SIZETYPE bwt_sz = 0;
    while (true)
    {
        fin.read((char *)SA, batch*sizeof(SIZETYPE));
        SIZETYPE gc = fin.gcount()/sizeof(SIZETYPE);
        if (gc==0)
            break;
        for (SIZETYPE i=0; i<gc; ++i, ++bwt_sz)
        {
            SIZETYPE j = (SA[i]>0 ? (SA[i]-1) : (revseq_sz-1));
            bwt[bwt_sz] = revseq[j];
        }
    }
    return bwt;
}

struct RankVec
{
    SIZETYPE** R1=NULL;
    R2TYPE** R2=NULL;
    const static SIZETYPE b1=8;
    const static SIZETYPE b2 = (SIZETYPE(std::numeric_limits<R2TYPE>::max()) + 1) / b1;
    const static NUCTYPE sigma=7;
    SIZETYPE C[sigma];
    NUCTYPE * bwt=NULL;
    SIZETYPE bwt_sz;

    RankVec()
    {
    }

    ~RankVec()
    {
        if (bwt)
            delete[] bwt;
        if(R1)
        {
            delete[] R1[1];
            delete[] R1;
        }
        if(R2)
        {
            delete[] R2[1];
            delete[] R2;
        }
    }

    template <typename U>
    void apply(U** & R, SIZETYPE block_size)
    {
        SIZETYPE sz=(bwt_sz-1)/block_size+2;
        if(R)
        {
            delete[] R[1];
            delete[] R;
        }
        R = new U*[sigma];
        R[1] = new U[(sigma-1)*sz];
        for(NUCTYPE c=2; c<sigma; ++c)
            R[c]=R[c-1]+sz;
    }

    void count()
    {
        apply(R1, b1*b2);
        apply(R2, b1);
        for(NUCTYPE c = 1; c < sigma; ++c)
        {
            R1[c][0] = 0;
            R2[c][0] = 0;
        }
        SIZETYPE j1 = 0;
        for(SIZETYPE i = 0, j2 = 0; i < bwt_sz; ++i)
        {
            if(i%(b1*b2) == 0)
            {
                ++j1;
                for(NUCTYPE c = 1; c < sigma; ++c)
                    R1[c][j1] = R1[c][j1-1];
            }
            if(i%b1==0)
            {
                ++j2;
                for(NUCTYPE c = 1; c < sigma; ++c)
                    if(j2%b2==0)
                        R2[c][j2]=0;
                    else
                        R2[c][j2]=R2[c][j2-1];
            }
            NUCTYPE c=bwt[i];
            if(c>0)
            {
                ++R1[c][j1];
                if(j2%b2!=0)
                    ++R2[c][j2];
            }
        }

        C[1]=1;
        for(NUCTYPE c = 2; c < sigma; ++c)
            C[c]=C[c-1]+R1[c-1][j1];
    }

    SIZETYPE rank(NUCTYPE c, SIZETYPE i)
    {
        SIZETYPE r=R1[c][i/(b1*b2)]+R2[c][i/b1];
        for(SIZETYPE j=(i/b1)*b1; j<i; ++j)
            if(bwt[j]==c)
                ++r;
        return r;
    }

    void saveRankVec(std::string bwtFile, std::string bwtRankFile)
    {
        std::ofstream fout(bwtFile);
        fout.write((char*)bwt, bwt_sz);
        fout.close();
        
        fout.open(bwtRankFile);
        fout.write((char*)C, sigma * sizeof(SIZETYPE));
        SIZETYPE sz=(bwt_sz+b1*b2-1)/(b1*b2)+1;
        fout.write((char*)R1[1], (sigma-1)*sz*sizeof(SIZETYPE));
        sz=(bwt_sz+b1-1)/b1+1;
        fout.write((char*)R2[1], (sigma-1)*sz*sizeof(R2TYPE));
    }

    void loadRankVec(std::string bwtFile, std::string bwtRankFile)
    {
        bwt_sz = std::filesystem::file_size(bwtFile);
        bwt = new NUCTYPE[bwt_sz];
        std::ifstream fin(bwtFile);
        fin.read((char*)bwt, bwt_sz);
        fin.close();

        apply(R1, b1*b2);
        apply(R2, b1);
        fin.open(bwtRankFile);
        fin.read((char*)C, sigma * sizeof(SIZETYPE));
        SIZETYPE sz=(bwt_sz+b1*b2-1)/(b1*b2)+1;
        fin.read((char*)R1[1], (sigma-1)*sz*sizeof(SIZETYPE));
        sz=(bwt_sz+b1-1)/b1+1;
        fin.read((char*)R2[1], (sigma-1)*sz*sizeof(R2TYPE));
    }
};


bool check_bwt(RankVec &bwtRank, NUCTYPE *revref)
{
    for(SIZETYPE i=0, j=bwtRank.bwt_sz-1; j>0; )
    {
        --j;
        if (bwtRank.bwt[i]!=revref[j])
            return false;
        NUCTYPE c=bwtRank.bwt[i];
        i=bwtRank.C[c]+bwtRank.rank(c,i);
    }
    return true;
}

#endif
