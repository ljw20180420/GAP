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
#include "Threadpool.h"

#define M64 1
#include "../external/gsacak.h"
#include "../external/gsacak.c"
#undef max

const static std::string int2base="#$NACGT";
static std::map<char, int8_t> base2int={{'#',0}, {'$',1}, {'N',2}, {'A',3}, {'C',4}, {'G',5}, {'T',6}, {'n',2}, {'a',3}, {'c',4}, {'g',5}, {'t',6}};

// std::tuple<uint8_t *, uint64_t, std::vector<std::pair<std::string, uint64_t>>> read_fasta_and_inverse(std::string REFfile)
// {
//     uint8_t *revref = (uint8_t *) malloc(std::filesystem::file_size(REFfile)*sizeof(uint8_t));
//     uint64_t revref_sz = 0;
//     std::vector<std::pair<std::string, uint64_t>> name_cumlen;
//     std::ifstream fin(REFfile);
//     for (std::string name, str; std::getline(std::getline(fin, name), str); )
//     {
//         name_cumlen.emplace_back(name, str.size());
//         for (uint64_t i=0; i<str.size(); ++i, ++revref_sz)
//             revref[revref_sz] = base2int[str[i]];
//         revref[revref_sz++] = 1;        
//     }
//     revref[revref_sz++] = 0;

//     // inverse the input sequence
//     std::reverse(revref, revref+revref_sz-2);
//     for(uint64_t i=0, j=name_cumlen.size()-1; i<j; ++i, --j)
//         std::swap(name_cumlen[i], name_cumlen[j]);

//     uint64_t tmp=name_cumlen[0].second;
//     name_cumlen[0].second=0;
//     for (uint64_t i=1; i<name_cumlen.size(); ++i)
//     {
//         std::swap(tmp, name_cumlen[i].second);
//         name_cumlen[i].second+=name_cumlen[i-1].second+1;
//     }
//     return std::make_tuple(revref, revref_sz, name_cumlen);
// }

void gsufsortSA(std::string SAfile, uint8_t *revref, uint64_t revref_sz)
{
    uint64_t *SA = (uint64_t *) malloc(revref_sz*sizeof(uint64_t));
    sacak(revref, SA, revref_sz);
    std::ofstream fSA(SAfile);
    fSA.write((char *)SA, revref_sz*sizeof(uint64_t));
    delete[] SA;
}

// void SA2GSA(std::string SAfile, std::string DAfile, std::string GAfile, uint64_t batch, std::vector<std::pair<std::string, uint64_t>> &name_cumlen)
// {
//     std::ifstream fin(SAfile);
//     std::ofstream fDA(DAfile), fGA(GAfile);
//     uint64_t *SA = (uint64_t *) malloc(batch*sizeof(uint64_t));
//     uint8_t *DA = (uint8_t *) malloc(batch*sizeof(uint8_t));
//     uint32_t *GA = (uint32_t *) malloc(batch*sizeof(uint32_t));
    
//     while (true)
//     {
//         fin.read((char *)SA, batch*sizeof(uint64_t));
//         uint64_t gc = fin.gcount()/sizeof(uint64_t);
//         if (gc==0)
//             break;
//         for (uint64_t i=0; i<gc; ++i)
//         {
//             uint64_t l=0, r=name_cumlen.size();
//             while (r-l>1)
//             {
//                 uint64_t m=(l+r)/2;
//                 if (name_cumlen[m].second>SA[i])
//                     r=m;
//                 else
//                     l=m;
//             }
//             DA[i] = l;
//             GA[i] = SA[i] - name_cumlen[l].second;
//         }
//         fDA.write((char *)DA, gc*sizeof(uint8_t));
//         fGA.write((char *)GA, gc*sizeof(uint32_t));
//     }
// }

uint8_t* SA2bwt(std::string SAfile, std::string BWTfile, uint64_t batch, uint8_t *revseq)
{
    uint64_t revseq_sz = std::filesystem::file_size(SAfile)/sizeof(uint64_t);
    std::ifstream fin(SAfile);
    std::ofstream fBWT(BWTfile);
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
    fBWT.write((char *)bwt, bwt_sz*sizeof(uint8_t));
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
    uint64_t vsz;

    RankVec(uint8_t *bwt, uint64_t vsz_) : vsz(vsz_)
    {
        apply_memory();
        std::cerr << "apply memory success\n"; // debug
        count(bwt);
        std::cerr << "count memory success\n"; // debug
    }

    RankVec(std::string bwtRankFile)
    {
        loadRankVec(bwtRankFile);
    }

    ~RankVec()
    {
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

    void apply_memory()
    {
        apply(R1, b1*b2);
        apply(R2, b1);
    }

    template <typename U>
    void apply(U** & R, uint64_t block_size)
    {
        uint64_t sz=(vsz-1)/block_size+2;
        release(R);
        R=new U*[sigma];
        R[0]=new U[sigma*sz];
        for(uint8_t c=1; c<sigma; ++c)
            R[c]=R[c-1]+sz;
    }

    void count(uint8_t *bwt)
    {
        for(uint8_t c=0; c<sigma; ++c)
        {
            R1[c][0]=0;
            R2[c][0]=0;
        }
        uint64_t j1=0;
        for(uint64_t i=0, j2=0; i<vsz; ++i)
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
        std::cerr << "R1 R2 success\n"; // debug
        std::cerr << j1 << '\n'; // debug

        C[0]=0;
        for(int c=1; c<sigma; ++c)
            C[c]=C[c-1]+R1[c-1][j1];
    }

    uint64_t rank(uint8_t c, uint64_t i, uint8_t *bwt)
    {
        if (i<0)
            return 0;
        uint64_t r=R1[c-1][i/(b1*b2)]+R2[c-1][i/b1];
        for(uint64_t j=(i/b1)*b1; j<=i; ++j)
            if(bwt[j]==c)
                ++r;
        return r;
    }

    void PreRange(int64_t& i, int64_t& j, int c, uint8_t *bwt)
    {
        i=C[c-1]+rank(c,i-1,bwt)+1;
        j=C[c-1]+rank(c,j,bwt);
    }

    void saveRankVec(std::string bwtRankFile)
    {
        std::ofstream fout(bwtRankFile);
        fout.write((char*)&vsz, sizeof(uint64_t));
        uint64_t sz=(vsz-1)/(b1*b2)+2;
        fout.write((char*)R1[0],this->sigma*sz*sizeof(R1[0][0]));
        sz=(vsz-1)/b1+2;
        fout.write((char*)R2[0],this->sigma*sz*sizeof(R2[0][0]));
    }

    void loadRankVec(std::string bwtRankFile)
    {
        std::ifstream fin(bwtRankFile);
        fin.read((char*)&vsz, sizeof(uint64_t));
        apply_memory();
        uint64_t sz=(vsz+b1*b2-1)/(b1*b2)+1;
        fin.read((char*)R1[0],this->sigma*sz*sizeof(R1[0][0]));
        sz=(vsz-1)/b1+2;
        fin.read((char*)R2[0],this->sigma*sz*sizeof(R2[0][0]));
    }
};


bool check_bwt(uint8_t *bwt, uint64_t bwt_sz, RankVec &bwtRank, uint8_t *revref)
{
    for(int64_t i=0, j=bwt_sz-2; j>=0; --j)
    {
        if (bwt[i]!=revref[j])
            return false;
        int c=bwt[i];
        i=bwtRank.C[c-1]+bwtRank.rank(c,i,bwt);
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

// struct BroWheel
// {
//     const static uint8_t sigma=6;
//     const static std::string int2base;
//     static std::map<char, int8_t> base2int;
//     bool reverse;
//     std::vector<uint8_t> sequence;
//     std::vector<uint8_t> bwt;
//     RankVec bwtRank = RankVec(sigma);
//     std::array<uint64_t,sigma> C;
//     const static uint8_t f=16;
//     std::vector<uint8_t> sra;
//     RankVec sraRank = RankVec(1);
//     uint64_t* ssa=NULL;
//     std::vector<std::pair<std::string, uint64_t>> name_cumlen;


//     ~BroWheel()
//     {
//         if(ssa)
//             delete[] ssa;
//     }

    // void readin(std::string file)
    // {
    //     sequence.clear();
    //     name_cumlen.clear();
    //     sequence.reserve(std::filesystem::file_size(file));
    //     std::ifstream fin(file);
    //     for (std::string name, str; std::getline(std::getline(fin, name), str); )
    //     {
    //         name_cumlen.emplace_back(name, str.size());
    //         for (uint64_t i=0; i<str.size(); ++i)
    //             sequence.push_back(base2int[str[i]]);
    //         sequence.push_back(1);            
    //     }
    //     sequence.push_back(0);

    //     // reverse the input sequence
    //     std::reverse(sequence.begin(), sequence.end()-2);
    //     for(uint64_t i=0, j=name_cumlen.size()-1; i<j; ++i, --j)
    //         std::swap(name_cumlen[i], name_cumlen[j]);

    //     uint64_t tmp=name_cumlen[0].second;
    //     name_cumlen[0].second=0;
    //     for (uint64_t i=1; i<name_cumlen.size(); ++i)
    //     {
    //         std::swap(tmp, name_cumlen[i].second);
    //         name_cumlen[i].second+=name_cumlen[i-1].second+1;
    //     }
    // }

    // void to_T(uint64_t h, uint64_t t, int64_t* SAI)
    // {
    //     uint8_t bsf=((sizeof(uint64_t)*CHAR_BIT-1)/3-1)*3;
    //     uint64_t i=t;
    //     SAI[i-h]=int64_t(-1)<<bsf;
    //     while(i>h)
    //     {
    //         --i;
    //         SAI[i-h]=(SAI[i-h+1]>>3)+(uint64_t(sequence[i])<<bsf);
    //     }
    // }

    // void index(int64_t ll, thread_pool & threads1, thread_pool & thread2)
    // {
    //     int64_t ssa_sz=(sequence.size()+f-1)/f;
    //     int64_t* SA=new int64_t[(threads1.size()+1)*(4*ll-1)+ssa_sz];
    //     int64_t* SAI=SA+2*ll-1;
    //     int64_t* SAmem=SAI+2*ll;
    //     int64_t* SAImem=SAmem+threads1.size()*(2*ll-1);
    //     int64_t* rssa=SAImem+threads1.size()*2*ll;
        
    //     bwt.resize(sequence.size());
    //     bwtRank.apply_memory(bwt.size());

    //     std::vector<int8_t> bwtt;
    //     bwtt.resize(sequence.size());
        
    //     int64_t N=sequence.size()/ll-1;
    //     std::queue<std::future<std::tuple<int64_t*,int64_t*,int64_t>>> futures1;
    //     int64_t m=(0>N-int64_t(threads1.size())?0:N-threads1.size());
    //     for(int64_t n=N-1; n>=m; --n)
    //        futures1.push(threads1.submit(std::bind(&BroWheel::SegmentSort, this, n, ll, SAmem, SAImem, threads1.size())));

    //     int64_t st=std::max(int64_t(0),N*ll);
    //     GetSA(st, sequence.size(), SA, SAI);
    //     bwt[SAI[0]] = 0;
    //     int64_t cutoff_size=sequence.size()-st;
    //     for(int64_t i=1; i<cutoff_size; ++i)
    //         bwt[SAI[i]] = sequence[st+i-1];
    //     int64_t R1_tail=bwtRank.count(cutoff_size, bwt);
    //     C[0]=0;
    //     for(int c=1; c<sigma; ++c)
    //         C[c]=C[c-1]+bwtRank.R1[c-1][R1_tail];

    //     for(int64_t n=N-1; n>=0; --n)
    //     {
    //         int64_t *SAn, *SAIn, pn;
    //         std::tie(SAn,SAIn,pn)=futures1.front().get();
    //         futures1.pop();
    //         std::future<void> future2=thread2.submit(std::bind(&BroWheel::RelativeOrder, this, SAI, SAIn, n, ll));

    //         for(int64_t i=0; i<pn; i+=2)
    //             SortSplit(SAn, SAIn[i], SAIn[i+1], NULL, SAI);
    //         future2.wait();
    //         bwt[SAI[0]] = sequence[(n+1)*ll-1];
    //         for(int64_t i=0; i<ll; ++i)
    //             SAI[SAn[i]]=SAIn[SAn[i]+ll]+i+1;
    //         cutoff_size+=ll;
    //         for(int64_t i=0, j=0, k=0; i<cutoff_size; ++i)
    //         {
    //             if(k>=ll || i!=SAI[SAn[k]])
    //                 bwtt[i] = bwt[j++];
    //             else
    //             {
    //                 if(SAn[k]>0)
    //                     bwtt[i] = sequence[n*ll+SAn[k]-1];
    //                 else
    //                     bwtt[i] = 0;
    //                 ++k;
    //             }
    //         }

    //         if(m>0)
    //         {
    //            --m;
    //            futures1.push(threads1.submit(std::bind(&BroWheel::SegmentSort, this, m, ll, SAmem, SAImem, threads1.size())));
    //         }
    //         std::swap(bwt,bwtt);
    //         R1_tail=bwtRank.count(cutoff_size, bwt);
    //         for(int c=1; c<sigma; ++c)
    //             C[c]=C[c-1]+bwtRank.R1[c-1][R1_tail];
    //     }
    //     bwtt.clear();

    //     sra.clear();
    //     sra.resize(sequence.size(), 0);
    //     sraRank.apply_memory(sra.size());
    //     for(int64_t k=ssa_sz-1, i=0, j=cutoff_size-2; k>=0; --k)
    //     {
    //         while(j>=k*f)
    //         {
    //             int c=bwt[i];
    //             i=C[c-1]+this->bwtRank.rank(c,i,bwt);
    //             --j;
    //         }
    //         rssa[k]=i;
    //         sra[i] = 1;
    //     }

    //     sraRank.count(sra.size(), sra);
    //     if(ssa)
    //         delete[] ssa;
    //     ssa=new int64_t[ssa_sz];
    //     for(int64_t k=0; k<ssa_sz; ++k)
    //         ssa[sraRank.rank(1, rssa[k], sra)-1]=k*f;
    //     delete[] SA;
    // }

    // bool check_bwt(int64_t cutoff_size)
    // {
    //     for(int64_t i=0, j=cutoff_size-2; j>=0; --j)
    //     {
    //         if (bwt[i]!=sequence[j])
    //             return false;
    //         int c=bwt[i];
    //         i=C[c-1]+this->bwtRank.rank(c,i,bwt);
    //     }
    //     return true;
    // }

    // void saveBroWheel(std::string file)
    // {
    //     std::ofstream fout(file+".idx", std::ios::binary);
    //     fout.write((char*)&reverse,sizeof(reverse));
    //     int64_t tmp=sequence.size();
    //     fout.write((char*)&tmp,sizeof(tmp));
    //     fout.write((char *)sequence.data(),sequence.size());

    //     int64_t name_cumlen_sz = name_cumlen.size();
    //     fout.write((char*)&name_cumlen_sz,sizeof(name_cumlen_sz));
    //     for (int i = 0; i < name_cumlen_sz; ++i)
    //     {
    //         int64_t str_sz = name_cumlen[i].first.size();
    //         fout.write((char*)&str_sz,sizeof(str_sz));
    //         fout.write(name_cumlen[i].first.data(),str_sz*sizeof(char));
    //         fout.write((char*)&name_cumlen[i].second,sizeof(name_cumlen[i].second));
    //     }

    //     tmp = bwt.size();
    //     fout.write((char*)&tmp,sizeof(tmp));
    //     fout.write((char*)bwt.data(),bwt.size());
    //     bwtRank.saveRankVec(fout, tmp);

    //     fout.write((char*)C.data(),sigma*sizeof(C[0]));

    //     tmp = sra.size();
    //     fout.write((char*)&tmp,sizeof(tmp));
    //     fout.write((char*)sra.data(),sra.size());
    //     sraRank.saveRankVec(fout, tmp);

    //     int64_t ssa_sz=(sequence.size()+f-1)/f;
    //     fout.write((char*)ssa,ssa_sz*sizeof(*ssa));
    //     fout.close();
    // }

    // void loadBroWheel(std::string file)
    // {
    //     if (!std::filesystem::exists(file+".idx"))
    //     {
    //         std::cerr << "long reference index " << file << ".idx does not exist";
    //         exit (EXIT_FAILURE);
    //     }

    //     std::ifstream fin(file+".idx", std::ios::binary);
    //     fin.read((char*)&reverse,sizeof(reverse));
    //     int64_t tmp;
    //     fin.read((char*)&tmp,sizeof(tmp));
    //     sequence.resize(tmp);
    //     fin.read((char *)sequence.data(),tmp);
    //     sequence.shrink_to_fit();

    //     int64_t name_cumlen_sz;
    //     fin.read((char*)&name_cumlen_sz,sizeof(name_cumlen_sz));
    //     name_cumlen.resize(name_cumlen_sz);
    //     for (int i = 0; i < name_cumlen_sz; ++i)
    //     {
    //         int64_t str_sz;
    //         fin.read((char*)&str_sz,sizeof(str_sz));
    //         name_cumlen[i].first.resize(str_sz);
    //         fin.read((char*)name_cumlen[i].first.data(),str_sz*sizeof(char));
    //         fin.read((char*)&name_cumlen[i].second,sizeof(name_cumlen[i].second));
    //     }

    //     fin.read((char*)&tmp,sizeof(tmp));
    //     bwt.resize(tmp);
    //     fin.read((char*)bwt.data(),tmp);
    //     bwtRank.loadRankVec(fin, tmp);

    //     fin.read((char*)C.data(),sigma*sizeof(C[0]));

    //     fin.read((char*)&tmp,sizeof(tmp));
    //     sra.resize(tmp);
    //     fin.read((char*)sra.data(),tmp);
    //     sraRank.loadRankVec(fin, tmp);

    //     int64_t ssa_sz=(sequence.size()+f-1)/f;
    //     if(ssa)
    //         delete[] ssa;
    //     ssa=new int64_t[ssa_sz];
    //     fin.read((char*)ssa,ssa_sz*sizeof(*ssa));
    //     fin.close();
    // }
    
    // void RelativeOrder(int64_t* SAI, int64_t* SAIn, int64_t n, int64_t ll)
    // {
    //     int c=sequence[(n+1)*ll-1];
    //     SAIn[2*ll-1]=C[c-1]+bwtRank.rank(c, SAI[0], bwt);
    //     for(int64_t i=ll-2; i>=0; --i)
    //     {
    //         c=sequence[n*ll+i];
    //         SAIn[i+ll]=C[c-1]+bwtRank.rank(c, SAIn[i+ll+1], bwt);
    //     }
    // }

    // std::tuple<int64_t*,int64_t*,int64_t> SegmentSort(int64_t n, int64_t ll, int64_t* SAmem, int64_t* SAImem, int para_sz)
    // {
    //     int64_t* SA=SAmem+(n%para_sz)*(2*ll-1);
    //     int64_t* SAI=SAImem+(n%para_sz)*2*ll;
    //     GetSA(n*ll, (n+2)*ll-1, SA, SAI);
    //     int64_t p=GetALS(n*ll, (n+2)*ll-1, SA, SAI, ll);
    //     return std::make_tuple(SA,SAI,p);
    // }
    
    // std::tuple<int64_t,int64_t*> OffSet(int* S, int64_t S_sz)
    // {
    //     int64_t i=0, j=bwt.size()-1;
    //     for(int64_t p=S_sz-1; p>=0; --p)
    //         PreRange(i,j,S[p]);
    //     int64_t O_sz=j-i+1;
    //     int64_t* O=new int64_t[O_sz];
    //     for(int64_t k=i; k<=j; ++k)
    //         O[k-i]=SimSuffix(k);
    //     return std::make_tuple(O_sz,O);
    // }
    
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
    
    // void PreRange(int64_t& i, int64_t& j, int c)
    // {
    //     i=C[c-1]+bwtRank.rank(c,i-1,bwt)+1;
    //     j=C[c-1]+bwtRank.rank(c,j,bwt);
    // }
    
    // int64_t GetALS(int64_t h, int64_t t, int64_t* SA, int64_t* SAI, int64_t ll)
    // {
    //     int64_t sz=t-h;
    //     for(int64_t i=0; i<ll; ++i)
    //         SA[SAI[i]]=i;
    //     for(int64_t i=0, j=0; i<sz; ++i)
    //         if(SA[i]>=0 && SA[i]<ll)
    //         {
    //             SA[j]=SA[i];
    //             SAI[SA[j]]=j;
    //             ++j;
    //         }
    //     for(int64_t i=0, k=0; i<ll; ++i)
    //         if(SAI[i]!=ll-1)
    //         {
    //             int64_t j=SA[SAI[i]+1];
    //             while(i+k<sz && j+k<sz && sequence[i+k+h]==sequence[j+k+h])
    //                 ++k;
    //             SAI[SAI[i]+ll]=k;
    //             if(k>0)
    //                 --k;
    //         }
    //     int64_t p=0;
    //     for(int64_t i=0; i<ll-1; ++i)
    //         if(SAI[i+ll]>=ll)
    //         {
    //             if(i==0 || SAI[i-1+ll]<ll)
    //                 SAI[p++]=i;
    //             if(i==ll-2 || SAI[i+1+ll]<ll)
    //                 SAI[p++]=i+1;
    //         }
    //     return p;
    // }

    // void GetSA(int64_t h, int64_t t, int64_t* SA, int64_t* SAI)
    // {
    //     int64_t sz=t-h;
    //     to_T(h,t,SAI);
    //     for(int64_t i=0; i<sz; ++i)
    //         SA[i]=i;
    //     SortSplit(SA, 0, sz-1, SAI, SAI);
    //     int64_t* SAIh=SAI+(sizeof(int64_t)*CHAR_BIT-1)/3;
    //     while(SA[0]>-sz)
    //     {
    //         int64_t i=0;
    //         int64_t k=0;
    //         do
    //         {
    //             if(SA[i]<0)
    //             {
    //                 k+=SA[i];
    //                 i-=SA[i];
    //             }
    //             else
    //             {
    //                 if(k<0)
    //                 {
    //                     SA[i+k]=k;
    //                     k=0;
    //                 }
    //                 int64_t j=SAI[SA[i]];
    //                 SortSplit(SA, i, j, SAI, SAIh);
    //                 i=j+1;
    //             }
    //         }
    //         while(i<sz);
    //         if(k<0)
    //             SA[i+k]=k;
    //         SAIh+=(SAIh-SAI);
    //     }
    // }

    // void SortSplit(int64_t* SA, int64_t p1, int64_t p2, int64_t* SAI, int64_t* SAIh)
    // {
    //     int64_t n=p2-p1+1;
    //     if(n<7)
    //     {
    //         InsertionSort(SA, p1, p2, SAI, SAIh);
    //         return;
    //     }
    //     int64_t m=p1+n/2;
    //     int64_t v=SAIh[SA[m]];
    //     if(n>7)
    //     {
    //         int64_t v1=SAIh[SA[p1]], v2=SAIh[SA[p2]];
    //         if(n>40)
    //         {
    //             int64_t s=n/8;
    //             v1=med3(v1,SAIh[SA[p1+s]],SAIh[SA[p1+2*s]]);
    //             v=med3(SAIh[SA[m-s]],v,SAIh[SA[m+s]]);
    //             v2=med3(SAIh[SA[p2-2*s]],SAIh[SA[p2-s]],v2);
    //         }
    //         v=med3(v1,v,v2);
    //     }
    //     int64_t pa=p1, pb=p1, pc=p2, pd=p2;
    //     while(true)
    //     {
    //         while(pb<=pc && SAIh[SA[pb]]<=v)
    //         {
    //             if(SAIh[SA[pb]]==v)
    //                 std::swap(SA[pa++],SA[pb]);
    //             ++pb;
    //         }
    //         while(pb<=pc && SAIh[SA[pc]]>=v)
    //         {
    //             if(SAIh[SA[pc]]==v)
    //                 std::swap(SA[pc],SA[pd--]);
    //             --pc;
    //         }
    //         if(pb>pc)
    //             break;
    //         std::swap(SA[pb++],SA[pc--]);
    //     }
    //     int64_t s=std::min(pa-p1,pb-pa);
    //     for(int64_t i=p1,j=pb-s; i<p1+s; ++i,++j)
    //         std::swap(SA[i],SA[j]);
    //     s=std::min(p2-pd,pd-pc);
    //     for(int64_t i=p2,j=pc+s; i>p2-s; --i,--j)
    //         std::swap(SA[i],SA[j]);
    //     pa=p1+pb-pa;
    //     pc=p2-pd+pc;
    //     if(pa>p1)
    //         SortSplit(SA, p1, pa-1, SAI, SAIh);
    //     if(SAI)
    //     {
    //         SAI[SA[pa]]=pc;
    //         if(pa==pc)
    //             SA[pa]=-1;
    //         else
    //             for (++pa; pa<=pc; ++pa)
    //                 SAI[SA[pa]]=pc;
    //     }
    //     if(pc<p2)
    //         SortSplit(SA, pc+1, p2, SAI, SAIh);
    // }

    // int64_t med3(int64_t a, int64_t b, int64_t c)
    // {
    //     if(a>b)
    //         if(b>c) return b;
    //         else if(a>c) return c;
    //             else return a;
    //     else
    //         if(b<c) return b;
    //         else if(a>c) return a;
    //         else return c;
    // }

    // void InsertionSort(int64_t* SA, int64_t p1, int64_t p2, int64_t* SAI, int64_t* SAIh)
    // {
    //     int64_t mv=std::numeric_limits<int64_t>::max();
    //     for (int64_t i=p1, j=i, k=j-1; i<=p2; mv=std::numeric_limits<int64_t>::max())
    //     {
    //         while (k<p2)
    //         {
    //             ++k;
    //             mv=std::min(mv, SAIh[SA[k]]);
    //         }

    //         while (j<=k)
    //             if (SAIh[SA[k]]==mv)
    //                 std::swap(SA[j++], SA[k]);
    //             else
    //                 --k;

    //         if (SAI)
    //             SAI[SA[i]]=j-1;

    //         if (j-i==1)
    //             SA[i++]=-1;
    //         else
    //             for (++i; i<j; ++i)
    //                 if (SAI)
    //                     SAI[SA[i]]=j-1;
    //     }
    // }

    // int64_t start_rev(int64_t start, int len)
    // {
    //     return sequence.size()-1-start-len;
    // }

    // std::tuple<std::string, int64_t> get_axis(int64_t s)
    // {
    //     s=sequence.size()-1-s;
    //     int64_t l=0, r=name_cumlen.size();
    //     while (r-l>1)
    //     {
    //         int64_t m=(l+r)/2;
    //         if (name_cumlen[m].second>s)
    //             r=m;
    //         else
    //             l=m;
    //     }
    //     s-=name_cumlen[l].second;
    //     if (r<name_cumlen.size())
    //         s=name_cumlen[r].second-name_cumlen[l].second-1-s;
    //     else
    //         s=sequence.size()-name_cumlen[l].second-1-s;
    //     return std::make_tuple(name_cumlen[l].first, s);
    // }
// };

// const std::string BroWheel::int2base="#$NACGT";
// std::map<char, int8_t> BroWheel::base2int={{'N',2}, {'A',3}, {'C',4}, {'G',5}, {'T',6}, {'n',2}, {'a',3}, {'c',4}, {'g',5}, {'t',6}};

#endif
