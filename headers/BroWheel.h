#ifndef BROWHEEL_H
#define BROWHEEL_H

#include <bits/stdc++.h>
#include <limits>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <tuple>
#include <list>
#include <map>
#include <cctype>
#include <algorithm>
#include <experimental/filesystem>
#include "Threadpool.h"

struct BinVec
{
    int sigma;
    int ww;
    int64_t sz;

    std::vector<char> BV;

    BinVec(int sigma_)
    {
        sigma=sigma_;
        sigma2ww();
    }

    void sigma2ww()
    {
        ww=0;
        int sigma_=sigma;
        while(sigma_>0)
        {
            sigma_/=2;
            ++ww;
        }
    }

    void resize(int64_t N)
    {
        BV.resize((ww*N+CHAR_BIT-1)/CHAR_BIT);
        sz=N;
    }
    
    bool getbit(int64_t j)
    {
        return BV[j/CHAR_BIT]<<(j%CHAR_BIT) & 128;
    }

    void setbit(int64_t j, bool logic)
    {
        if(logic)
            BV[j/CHAR_BIT]|=128>>(j%CHAR_BIT);
        else
            BV[j/CHAR_BIT]&=~(128>>(j%CHAR_BIT));
    }

    int64_t operator()(int64_t i)
    {
        int64_t j=ww*i;
        int64_t c=getbit(j);
        while(j<ww*(i+1)-1)
        {
            c<<=1;
            c+=getbit(++j);
        }
        return c;
    }
    
    void set(int64_t i, int c)
    {
        int64_t j=ww*(i+1);
        while(j>ww*i)
        {
            setbit(--j, c & 1);
            c/=2;
        }
    }

    void push_back(int c)
    {
        resize(++sz);
        set(sz-1,c);
    }

    void clear()
    {
        BV.clear();
        sz=0;
    }

    void reserve(int64_t N)
    {
        BV.reserve((N*ww+CHAR_BIT-1)/CHAR_BIT);
    }

    int64_t size()
    {
        return sz;
    }

    void shrink_to_fit()
    {
        BV.resize((sz*ww+CHAR_BIT-1)/CHAR_BIT);
        BV.shrink_to_fit();
    }

    void saveBinVec(std::ofstream & fout)
    {
        fout.write((char*)&sigma,sizeof(sigma));
        fout.write((char*)&sz,sizeof(sz));
        int64_t tmp=BV.size();
        fout.write((char*)&tmp,sizeof(tmp));
        fout.write(BV.data(),BV.size());
    }

    void loadBinVec(std::ifstream & fin)
    {
        fin.read((char*)&sigma,sizeof(sigma));
        sigma2ww();
        fin.read((char*)&sz,sizeof(sz));
        int64_t tmp;
        fin.read((char*)&tmp,sizeof(tmp));
        BV.resize(tmp);
        fin.read(BV.data(),tmp);
        BV.shrink_to_fit();
    }
};

struct RankVec : BinVec
{
    int64_t** R1=NULL;
    uint16_t** R2=NULL;
    const static int b1=8;
    const static int b2=256*32;

    RankVec(int sigma_) : BinVec(sigma_)
    {}

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
    void apply(U** & R, int block_size)
    {
        int64_t sz=(this->size()-1)/block_size+2;
        release(R);
        R=new U*[this->sigma];
        R[0]=new U[this->sigma*sz];
        for(int i=1; i<this->sigma; ++i)
            R[i]=R[i-1]+sz;
    }

    int64_t count(int64_t cutoff_size)
    {
        for(int c=0; c<this->sigma; ++c)
        {
            R1[c][0]=0;
            R2[c][0]=0;
        }
        int64_t j1=0;
        for(int64_t i=0, j2=0; i<cutoff_size; ++i)
        {
            if(i%(b1*b2)==0)
            {
                ++j1;
                for(int c=0; c<this->sigma; ++c)
                    R1[c][j1]=R1[c][j1-1];
            }
            if(i%b1==0)
            {
                ++j2;
                for(int c=0; c<this->sigma; ++c)
                    if(j2%b2==0)
                        R2[c][j2]=0;
                    else
                        R2[c][j2]=R2[c][j2-1];
            }
            int c=(*this)(i);
            if(c>0)
            {
                ++R1[c-1][j1];
                if(j2%b2!=0)
                    ++R2[c-1][j2];
            }
        }
        return j1;
    }

    int64_t rank(int c, int64_t i)
    {
        if (i<0)
            return 0;
        int64_t r=R1[c-1][i/(b1*b2)]+R2[c-1][i/b1];
        for(int64_t j=(i/b1)*b1; j<=i; ++j)
            if((*this)(j)==c)
                ++r;
        return r;
    }

    void saveRankVec(std::ofstream & fout)
    {
        this->saveBinVec(fout);
        int64_t sz=(this->size()-1)/(b1*b2)+2;
        fout.write((char*)R1[0],this->sigma*sz*sizeof(R1[0][0]));
        sz=(this->size()-1)/b1+2;
        fout.write((char*)R2[0],this->sigma*sz*sizeof(R2[0][0]));
    }

    void loadRankVec(std::ifstream & fin)
    {
        this->loadBinVec(fin);
        apply_memory();
        int64_t sz=(this->size()-1)/(b1*b2)+2;
        fin.read((char*)R1[0],this->sigma*sz*sizeof(R1[0][0]));
        sz=(this->size()-1)/b1+2;
        fin.read((char*)R2[0],this->sigma*sz*sizeof(R2[0][0]));
    }
};

struct BroWheel
{
    const static int block_size=64*1024;
    const static int sigma=6;
    const static std::string int2base;
    BinVec sequence=BinVec(sigma);
    RankVec bwt=RankVec(sigma);
    std::array<int64_t,sigma> C;
    const static int f=16;
    RankVec sra=RankVec(1);
    int64_t* ssa=NULL;

    ~BroWheel()
    {
        if(ssa)
            delete[] ssa;
    }

    void readin(std::list<std::string> files, bool reverse)
    {
        int64_t sequence_sz=0;
        std::experimental::filesystem::path path=std::experimental::filesystem::current_path();
        for(auto file : files)
            sequence_sz+=std::experimental::filesystem::file_size(path/file);
        sequence.clear();
        sequence.reserve(sequence_sz);
        
        char* str=new char[block_size];
        for(auto file : files)
        {
            bool flag=false;
            std::ifstream fin(file, std::ifstream::binary);
            do
            {
                fin.read(str,block_size);
                for(int i=0; i<(fin.good()?block_size:fin.gcount()); ++i)
                {
                    if(!isspace(str[i]))
                    {
                        if(!flag)
                            flag=true;
                        sequence.push_back((str[i]+6)>>1);
                    }
                    else
                    {
                        if(flag)
                        {
                            sequence.push_back(1);
                            flag=false;
                        }
                    }
                }
            }while(fin.good());
            if(flag)
            {
                sequence.push_back(1);
            }
            fin.close();
        }
        sequence.set(sequence.size()-1,0);
        sequence.shrink_to_fit();
        if(reverse)
            for(int64_t i=0,j=sequence.size()-2; i<j; ++i,--j)
            {
                int c=sequence(i);
                sequence.set(i,sequence(j));
                sequence.set(j,c);
            }
        delete[] str;
    }

    void to_T(int64_t h, int64_t t, int64_t* SAI)
    {
        int64_t bsf=((sizeof(int64_t)*CHAR_BIT-1)/sequence.ww-1)*sequence.ww;
        int64_t i=t-1;
        SAI[i-h]=sequence(i)<<bsf;
        while(i>h)
        {
            --i;
            SAI[i-h]=(SAI[i-h+1]>>sequence.ww)+(sequence(i)<<bsf);
        }
    }

    void index(int64_t ll, thread_pool & threads1, thread_pool & thread2)
    {
        int64_t ssa_sz=(sequence.size()+f-1)/f;
        int64_t* SA=new int64_t[(threads1.size()+1)*(4*ll-1)+ssa_sz];
        int64_t* SAI=SA+2*ll-1;
        int64_t* SAmem=SAI+2*ll;
        int64_t* SAImem=SAmem+threads1.size()*(2*ll-1);
        int64_t* rssa=SAImem+threads1.size()*2*ll;
        if(ssa)
            delete[] ssa;
        ssa=new int64_t[ssa_sz];
        bwt.resize(sequence.size());
        bwt.apply_memory();
        sra.resize(sequence.size());
        sra.apply_memory();
        BinVec bwtt(sigma);
        bwtt.resize(sequence.size());
        
        int64_t N=sequence.size()/ll-1;
        std::queue<std::future<std::tuple<int64_t*,int64_t*,int64_t>>> futures1;
        int64_t m=(0>N-int64_t(threads1.size())?0:N-threads1.size());
        for(int64_t n=N-1; n>=m; --n)
           futures1.push(threads1.submit(std::bind(&BroWheel::SegmentSort, this, n, ll, SAmem, SAImem, threads1.size())));
        
        int64_t st=std::max(int64_t(0),N*ll);
        GetSA(st, sequence.size(), SA, SAI);
        bwt.set(SAI[0],0);
        int64_t cutoff_size=sequence.size()-st;
        for(int64_t i=1; i<cutoff_size; ++i)
            bwt.set(SAI[i],sequence(st+i-1));
        int64_t R1_tail=bwt.count(cutoff_size);
        C[0]=0;
        for(int c=1; c<sigma; ++c)
            C[c]=C[c-1]+bwt.R1[c-1][R1_tail];
        
        for(int64_t n=N-1; n>=0; --n)
        {
            int64_t *SAn, *SAIn, pn;
            std::tie(SAn,SAIn,pn)=futures1.front().get();
            futures1.pop();
            std::future<void> future2=thread2.submit(std::bind(&BroWheel::RelativeOrder, this, SAI, SAIn, n, ll));

            for(int64_t i=0; i<pn; i+=2)
                SortSplit(SAn, SAIn[i], SAIn[i+1], NULL, SAI);
            future2.wait();
            bwt.set(SAI[0],sequence((n+1)*ll-1));
            for(int64_t i=0; i<ll; ++i)
                SAI[SAn[i]]=SAIn[SAn[i]+ll]+i+1;
            cutoff_size+=ll;
            for(int64_t i=0, j=0, k=0; i<cutoff_size; ++i)
            {
                if(k>=ll || i!=SAI[SAn[k]])
                    bwtt.set(i,bwt(j++));
                else
                {
                    if(SAn[k]>0)
                        bwtt.set(i,sequence(n*ll+SAn[k]-1));
                    else
                        bwtt.set(i,0);
                    ++k;
                }
            }
            if(m>0)
            {
               --m;
               futures1.push(threads1.submit(std::bind(&BroWheel::SegmentSort, this, m, ll, SAmem, SAImem, threads1.size())));
            }
            std::swap(bwt.BV,bwtt.BV);
            R1_tail=bwt.count(cutoff_size);
            for(int c=1; c<sigma; ++c)
                C[c]=C[c-1]+bwt.R1[c-1][R1_tail];
        }

        for(int64_t k=ssa_sz-1, i=0, j=cutoff_size-2; k>=0; --k)
        {
            while(j>=k*f)
            {
                int c=bwt(i);
                i=C[c-1]+this->bwt.rank(c,i);
                --j;
            }
            rssa[k]=i;
            sra.set(i,1);
        }
        sra.count(sra.size());
        for(int64_t k=0; k<ssa_sz; ++k)
            ssa[sra.rank(1,rssa[k])-1]=k*f;
        delete[] SA;
    }

    void saveBroWheel(std::string file)
    {
        std::ofstream fout(file, std::ios::binary);
        sequence.saveBinVec(fout);
        bwt.saveRankVec(fout);
        fout.write((char*)C.data(),sigma*sizeof(C[0]));
        sra.saveRankVec(fout);
        int64_t ssa_sz=(sequence.size()+f-1)/f;
        fout.write((char*)ssa,ssa_sz*sizeof(*ssa));
        fout.close();
    }

    void loadBroWheel(std::string file)
    {
        std::ifstream fin(file, std::ios::binary);
        sequence.loadBinVec(fin);
        bwt.loadRankVec(fin);
        fin.read((char*)C.data(),sigma*sizeof(C[0]));
        sra.loadRankVec(fin);
        int64_t ssa_sz=(sequence.size()+f-1)/f;
        if(ssa)
            delete[] ssa;
        ssa=new int64_t[ssa_sz];
        fin.read((char*)ssa,ssa_sz*sizeof(*ssa));
        fin.close();
    }
    
    void RelativeOrder(int64_t* SAI, int64_t* SAIn, int64_t n, int64_t ll)
    {
        int c=sequence((n+1)*ll-1);
        SAIn[2*ll-1]=C[c-1]+bwt.rank(c,SAI[0]);
        for(int64_t i=ll-2; i>=0; --i)
        {
            c=sequence(n*ll+i);
            SAIn[i+ll]=C[c-1]+bwt.rank(c,SAIn[i+ll+1]);
        }
    }

    std::tuple<int64_t*,int64_t*,int64_t> SegmentSort(int64_t n, int64_t ll, int64_t* SAmem, int64_t* SAImem, int para_sz)
    {
        int64_t* SA=SAmem+(n%para_sz)*(2*ll-1);
        int64_t* SAI=SAImem+(n%para_sz)*2*ll;
        GetSA(n*ll, (n+2)*ll-1, SA, SAI);
        int64_t p=GetALS(n*ll, (n+2)*ll-1, SA, SAI, ll);
        return std::make_tuple(SA,SAI,p);
    }
    
    std::tuple<int64_t,int64_t*> OffSet(int* S, int64_t S_sz)
    {
        int64_t i=0, j=bwt.size()-1;
        for(int64_t p=S_sz-1; p>=0; --p)
            PreRange(i,j,S[p]);
        int64_t O_sz=j-i+1;
        int64_t* O=new int64_t[O_sz];
        for(int64_t k=i; k<=j; ++k)
            O[k-i]=SimSuffix(k);
        return std::make_tuple(O_sz,O);
    }
    
    int64_t SimSuffix(int64_t i)
    {
        int j=0;
        while(sra(i)==0)
        {
            int c=bwt(i);
            i=C[c-1]+bwt.rank(c,i);
            ++j;
        }
        return ssa[sra.rank(1,i)-1]+j;
    }
    
    void PreRange(int64_t& i, int64_t& j, int c)
    {
        i=C[c-1]+bwt.rank(c,i-1)+1;
        j=C[c-1]+bwt.rank(c,j);
    }
    
    int64_t GetALS(int64_t h, int64_t t, int64_t* SA, int64_t* SAI, int64_t ll)
    {
        int64_t sz=t-h;
        for(int64_t i=0; i<ll; ++i)
            SA[SAI[i]]=i;
        for(int64_t i=0, j=0; i<sz; ++i)
            if(SA[i]>=0 && SA[i]<ll)
            {
                SA[j]=SA[i];
                SAI[SA[j]]=j;
                ++j;
            }
        for(int64_t i=0, k=0; i<ll; ++i)
            if(SAI[i]!=ll-1)
            {
                int64_t j=SA[SAI[i]+1];
                while(i+k<sz && j+k<sz && sequence(i+k+h)==sequence(j+k+h))
                    ++k;
                SAI[SAI[i]+ll]=k;
                if(k>0)
                    --k;
            }
        int64_t p=0;
        for(int64_t i=0; i<ll-1; ++i)
            if(SAI[i+ll]>=ll)
            {
                if(i==0 || SAI[i-1+ll]<ll)
                    SAI[p++]=i;
                if(i==ll-2 || SAI[i+1+ll]<ll)
                    SAI[p++]=i+1;
            }
        return p;
    }

    void GetSA(int64_t h, int64_t t, int64_t* SA, int64_t* SAI)
    {
        int64_t sz=t-h;
        to_T(h,t,SAI);
        for(int64_t i=0; i<sz; ++i)
            SA[i]=i;
        SAI[sz]=0;
        SortSplit(SA, 0, sz-1, SAI, SAI);
        int64_t* SAIh=SAI+(sizeof(int64_t)*CHAR_BIT-1)/sequence.ww;
        while(SA[0]>-sz)
        {
            int64_t i=0;
            int64_t k=0;
            do
            {
                if(SA[i]<0)
                {
                    k+=SA[i];
                    i-=SA[i];
                }
                else
                {
                    if(k<0)
                    {
                        SA[i+k]=k;
                        k=0;
                    }
                    int64_t j=SAI[SA[i]];
                    SortSplit(SA, i, j, SAI, SAIh);
                    i=j+1;
                }
            }
            while(i<sz);
            if(k<0)
                SA[i+k]=k;
            SAIh+=(SAIh-SAI);
        }
    }

    void SortSplit(int64_t* SA, int64_t p1, int64_t p2, int64_t* SAI, int64_t* SAIh)
    {
        int64_t n=p2-p1+1;
        if(n<7)
        {
            InsertionSort(SA, p1, p2, SAI, SAIh);
            return;
        }
        int64_t m=p1+n/2;
        int64_t v=SAIh[SA[m]];
        if(n>7)
        {
            int64_t v1=SAIh[SA[p1]], v2=SAIh[SA[p2]];
            if(n>40)
            {
                int64_t s=n/8;
                v1=med3(v1,SAIh[SA[p1+s]],SAIh[SA[p1+2*s]]);
                v=med3(SAIh[SA[m-s]],v,SAIh[SA[m+s]]);
                v2=med3(SAIh[SA[p2-2*s]],SAIh[SA[p2-s]],v2);
            }
            v=med3(v1,v,v2);
        }
        int64_t pa=p1, pb=p1, pc=p2, pd=p2;
        while(true)
        {
            while(pb<=pc && SAIh[SA[pb]]<=v)
            {
                if(SAIh[SA[pb]]==v)
                    std::swap(SA[pa++],SA[pb]);
                ++pb;
            }
            while(pb<=pc && SAIh[SA[pc]]>=v)
            {
                if(SAIh[SA[pc]]==v)
                    std::swap(SA[pc],SA[pd--]);
                --pc;
            }
            if(pb>pc)
                break;
            std::swap(SA[pb++],SA[pc--]);
        }
        int64_t s=std::min(pa-p1,pb-pa);
        for(int64_t i=p1,j=pb-s; i<p1+s; ++i,++j)
            std::swap(SA[i],SA[j]);
        s=std::min(p2-pd,pd-pc);
        for(int64_t i=p2,j=pc+s; i>p2-s; --i,--j)
            std::swap(SA[i],SA[j]);
        pa=p1+pb-pa;
        pc=p2-pd+pc;
        if(pa>p1)
            SortSplit(SA, p1, pa-1, SAI, SAIh);
        if(SAI)
        {
            SAI[SA[pa]]=pc;
            if(pa==pc)
                SA[pa]=-1;
            else
                do
                    SAI[SA[++pa]]=pc;
                while(pa!=pc);
        }
        if(pc<p2)
            SortSplit(SA, pc+1, p2, SAI, SAIh);
    }

    int64_t med3(int64_t a, int64_t b, int64_t c)
    {
        if(a>b)
            if(b>c) return b;
            else if(a>c) return c;
                else return a;
        else
            if(b<c) return b;
            else if(a>c) return a;
            else return c;
    }

    void InsertionSort(int64_t* SA, int64_t p1, int64_t p2, int64_t* SAI, int64_t* SAIh)
    {
        for(int64_t i=p1+1; i<=p2; ++i)
        {
            int64_t j=i;
            while(j>p1 && SAIh[SA[j]]<SAIh[SA[j-1]])
            {
                std::swap(SA[j],SA[j-1]);
                --j;
            }
        }
        if(SAI)
        {
            int64_t k=p1;
            int64_t v=SAIh[SA[k]];
            for(int64_t i=p1+1; i<=p2; ++i)
            {
                int64_t u=SAIh[SA[i]];
                if(u>v)
                {
                    for(int64_t j=k; j<i; ++j)
                        SAI[SA[j]]=i-1;
                    if(i-k==1)
                        SA[k]=-1;
                    k=i;
                    v=u;
                }
            }
            for(int64_t j=k; j<=p2; ++j)
                SAI[SA[j]]=p2;
            if(k==p2)
                SA[k]=-1;
        }
    }

    int64_t start_rev(int64_t start, int len)
    {
        return sequence.size()-1-start-len;
    }

    static int base2int(char c)
    {
        return (c + 6) >> 1 & 7;
    }
};

std::string readin_local(std::list<std::string> files, bool reverse)
{
    std::string str, tmp;
    for (auto &file : files)
    {
        std::ifstream fin(file);
        fin >> tmp;
        while (fin.good())
        {
            str+=tmp;
            str.push_back('L');
            fin >> tmp;
        }
        if (!empty(str))
        {
            str.pop_back();
            str.push_back('K');
        }
    }
    if (reverse)
        for (size_t i=0, j=str.size()-2; i<j; ++i, --j)
            std::swap(str[i],str[j]);
    return str;
}

const std::string BroWheel::int2base="KLNACTG";

#endif
