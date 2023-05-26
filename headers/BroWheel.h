#ifndef BROWHEEL_H
#define BROWHEEL_H

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

struct Memory
{
    size_t chunk_sz;
    std::vector<char*> chunks;
    char* now;
    size_t remain;
    size_t chunk_id;
    
    ~Memory()
    {
        for(size_t i=0; i<chunks.size(); ++i)
            delete[] chunks[i];
    }
    
    void Initial(size_t chunk_sz_)
    {
        for(size_t i=0; i<chunks.size(); ++i)
            delete[] chunks[i];
        chunks.clear();
        chunk_sz=chunk_sz_;
        chunks.emplace_back(new char[chunk_sz]);
        now=chunks.back();
        remain=chunk_sz;
        chunk_id=0;
    }
    
    template <typename U>
    U* heap_alloc(size_t num)
    {
        size_t size=num*sizeof(U);
        void* ptr=now;
        if(!(now=(char*)std::align(alignof(U),size,ptr,remain)))
        {
            ++chunk_id;
            if(chunk_id==chunks.size())
            {
                chunks.emplace_back(new char[chunk_sz]);
                now=chunks.back();
            }
            else
                now=chunks[chunk_id];
            remain=chunk_sz;
            ptr=now;
            now=(char*)std::align(alignof(U),size,ptr,remain);
        }
        now+=size;
        remain-=size;
        return (U*)(now-size);
    }
    
    void clear()
    {
        chunk_id=0;
        now=chunks.front();
        remain=chunk_sz;
    }
};

template <typename T>
struct BinWord
{
    static const int ww=3;
    static const int r=(sizeof(T)*8-2)/ww;
    static const int bsf=(r-1)*ww;
    static const int block_size=64*1024;

    std::vector<bool> BV;
    
    void readin(std::list<std::string> & files, bool reverse)
    {
        T G_sz=0;
        std::experimental::filesystem::path path=std::experimental::filesystem::current_path();
        for(auto file : files)
            G_sz+=std::experimental::filesystem::file_size(path/file);
        BV.clear();
        BV.reserve(G_sz*ww);
        
        char* str=new char[block_size];
        for(auto file : files)
        {
            bool flag=false;
            std::ifstream fin(file, std::ifstream::binary);
            do
            {
                fin.read(str,block_size);
                for(int i=0, tmp=fin.good()?block_size:fin.gcount(); i<tmp; ++i)
                {
                    if(!isspace(str[i]))
                    {
                        if(!flag)
                            flag=true;
                        int tmp=str[i]+6;
                        BV.push_back(tmp & 8);
                        BV.push_back(tmp & 4);
                        BV.push_back(tmp & 2);
                    }
                    else
                    {
                        if(flag)
                        {
                            BV.push_back(0);
                            BV.push_back(0);
                            BV.push_back(1);
                            flag=false;
                        }
                    }
                }
            }while(fin.good());
            if(flag)
            {
                BV.push_back(0);
                BV.push_back(0);
                BV.push_back(1);
            }
            fin.close();
        }
        BV.back()=0;
        BV.shrink_to_fit();
        if(reverse)
            for(T i=0,j=BV.size()-6; i<j; j-=6)
            {
                std::swap(BV[i++],BV[j++]);
                std::swap(BV[i++],BV[j++]);
                std::swap(BV[i++],BV[j++]);
            }
        delete[] str;
    }
    
    void resize(T N)
    {
        BV.resize(ww*N);
        BV.shrink_to_fit();
    }
    
    T operator()(T i)
    {
        T j=ww*i, v=BV[j]<<1;
        if(BV[++j]) ++v;
        v<<=1;
        if(BV[++j]) ++v;
        return v;
    }
    
    void set(T i, T v)
    {
        T j=ww*i;
        BV[j]=(v & 4);
        BV[++j]=(v & 2);
        BV[++j]=(v & 1);
    }
    
    void to_T(T h, T t, T* SAI)
    {
        T i=t-1, k=i-h;
        SAI[k]=(*this)(i)<<bsf;
        while(k>0)
        {
            --k;
            SAI[k]=(SAI[k+1]>>ww)+((*this)(--i)<<bsf);
        }
    }
    
    T size()
    {
        return BV.size()/ww;
    }
    
    void to_chars(char* bytes)
    {
        T sz=(BV.size()+sizeof(char)-1)/sizeof(char);
        T j=0;
        for(T i=0; i<sz; ++i)
        {
            bytes[i]=0;
            for(T k=sizeof(char)-1; (k>-1)&&(j<BV.size()); --k,++j)
                bytes[i]|= BV[j]<<k;
        }
    }
    
    void from_chars(char* bytes, T BV_sz)
    {
        BV.resize(BV_sz);
        BV.shrink_to_fit();
        T sz=(BV.size()+sizeof(char)-1)/sizeof(char);
        T j=0;
        for(T i=0; i<sz; ++i)
            for(T k=7; (k>-1)&&(j<BV.size()); --k,++j)
                BV[j]=bytes[i] & 1UL<<k;
    }
};


template <typename T>
struct BroWheel : BinWord<T>, Memory
{
    static const int para_sz=6;
    static const int sigma=6;
    static const int oo=256;
    static const int oo2=oo*oo;
    static const int f=32;
    
    T* C;
    T** Occ1;
    uint16_t** Occ2;
    T* SSA;
    uint16_t** RSRA;
    
    void index(BinWord<T> & G, T ll, thread_pool & threads1, thread_pool & thread2, std::string outfile)
    {
        T RSRA_sz=(G.size()+oo2-1)/oo2+1, SSA_sz=(G.size()+f-2)/f, Occ1_sz=(G.size()+oo2-1)/oo2+1, Occ2_sz=(G.size()+oo-1)/oo+1;
        T* SA=new T[(para_sz+1)*(4*ll-1)+SSA_sz];
        T* SAI=SA+2*ll-1;
        T* SAmem=SAI+2*ll;
        T* SAImem=SAmem+para_sz*(2*ll-1);
        T* SRA=SAImem+para_sz*2*ll;
        T** ptr=new T*[RSRA_sz];
        Initial((sigma*(Occ1_sz+1)+SSA_sz+3)*sizeof(T)+(sigma+1)*sizeof(T*)+(sigma*Occ2_sz+SSA_sz+2)*sizeof(uint16_t)+(sigma+RSRA_sz+2)*sizeof(uint16_t*));
        C=heap_alloc<T>(sigma);
        Occ1=heap_alloc<T*>(sigma);
        Occ2=heap_alloc<uint16_t*>(sigma);
        RSRA=heap_alloc<uint16_t*>(RSRA_sz);
        RSRA[0]=heap_alloc<uint16_t>(SSA_sz);
        SSA=heap_alloc<T>(SSA_sz);
        Occ1[0]=heap_alloc<T>(sigma*Occ1_sz);
        Occ2[0]=heap_alloc<uint16_t>(sigma*Occ2_sz);
        for(T i=1; i<sigma; ++i)
        {
            Occ1[i]=&Occ1[i-1][Occ1_sz];
            Occ2[i]=&Occ2[i-1][Occ2_sz];
        }
        this->resize(G.size());
        BinWord<T> BWTT;
        BWTT.resize(G.size());
        
        T N=G.size()/ll-1;
        std::queue<std::future<std::tuple<T*,T*,T>>> futures1;
        T m=std::max(T(0),N-para_sz);
        for(T n=N-1; n>=m; --n)
            futures1.push(threads1.submit(std::bind(&BroWheel::SegmentSort, this, G, n, ll, SAmem, SAImem)));
        
        int st=std::max(T(0),N*ll);
        GetSA(G, st, G.size(), SA, SAI);
        this->set(SAI[0],0);
        T BWT_sz=G.size()-st;
        for(T i=1; i<BWT_sz; ++i)
            this->set(SAI[i],G(st+i-1));
        for(T i=0; i<sigma; ++i)
        {
            C[i]=0;
            Occ1[i][0]=0;
            Occ2[i][0]=0;
        }
        CountBWT(BWT_sz);
        
        for(T n=N-1; n>=0; --n)
        {
            T *SAn, *SAIn, pn;
            std::tie(SAn,SAIn,pn)=futures1.front().get();
            futures1.pop();
            std::future<void> future2=thread2.submit(std::bind(&BroWheel::RelativeOrder, this, G, SAI, SAIn, n, ll));
            for(T i=0; i<pn; i+=2)
                SortSplit(SAn, SAIn[i], SAIn[i+1], SAI, SAI, false);
            future2.wait();
            this->set(SAI[0],G((n+1)*ll-1));
            for(T i=0; i<ll; ++i)
                SAI[SAn[i]]=SAIn[SAn[i]+ll]+i+1;
            BWT_sz+=ll;
            for(T i=0, j=0, k=0; i<BWT_sz; ++i)
            {
                if(k>=ll || i!=SAI[SAn[k]])
                    BWTT.set(i,(*this)(j++));
                else
                {
                    if(SAn[k]>0)
                        BWTT.set(i,G(n*ll+SAn[k]-1));
                    else
                        BWTT.set(i,0);
                    ++k;
                }
            }
            if(m>0)
            {
                --m;
                futures1.push(threads1.submit(std::bind(&BroWheel::SegmentSort, this, G, m, ll, SAmem, SAImem)));
            }
            std::swap(this->BV,BWTT.BV);
            CountBWT(BWT_sz);
        }
        
        for(T p=1; p<RSRA_sz-1; ++p)
            RSRA[p]=RSRA[0];
        T i=0;
        for(T j=BWT_sz-2, tmp; j>=(SSA_sz-1)*f; --j)
        {
            tmp=(*this)(i);
            i=C[tmp-1]+this->GetOcc(i,tmp);
        }
        SRA[SSA_sz-1]=i;
        for(T k=SSA_sz-2; k>=0; --k)
        {
            for(T c=0, tmp; c<f; ++c)
            {
                tmp=(*this)(i);
                i=C[tmp-1]+this->GetOcc(i,tmp);
            }
            SRA[k]=i;
            ++RSRA[i/oo2+1];
        }
        ptr[0]=ptr[1]=&SSA[0];
        for(T p=1, tmp; p<RSRA_sz-1; ++p)
        {
            tmp=RSRA[p]-RSRA[0];
            ptr[p+1]=ptr[p]+tmp;
            RSRA[p]=RSRA[p-1]+tmp;
        }
        RSRA[RSRA_sz-1]=&RSRA[0][SSA_sz];
        for(T k=0, p; k<SSA_sz; ++k)
        {
            p=SRA[k]/oo2+1;
            *(ptr[p]++)=k;
        }
        for(T p=1, *pptr=ptr[0]; p<RSRA_sz; pptr=ptr[p],++p)
        {
            std::sort(pptr, ptr[p],[SRA](T i1, T i2) {return SRA[i1]<SRA[i2];});
        }
        for(T k=0; k<SSA_sz; ++k)
        {
            RSRA[0][k]=SRA[SSA[k]]%oo2;
            SSA[k]*=f;
        }
        delete[] SA;
        delete[] ptr;
        
        std::ofstream fout(outfile, std::ios::binary);
        T BV_sz=G.BV.size();
        fout.write((char*)&BV_sz,sizeof(T));
        T sz=(BV_sz+7)/8;
        char* bytes=new char[sz];
        G.to_chars(bytes);
        fout.write(bytes,sz);
        this->to_chars(bytes);
        fout.write(bytes,sz);
        delete[] bytes;
        fout.write((char*)&RSRA_sz,sizeof(T));
        fout.write((char*)&SSA_sz,sizeof(T));
        fout.write((char*)&Occ1_sz,sizeof(T));
        fout.write((char*)&Occ2_sz,sizeof(T));
        for(T i=1; i<RSRA_sz; ++i)
        {
            T tmp=RSRA[i]-RSRA[0];
            fout.write((char*)&tmp,sizeof(T));
        }
        fout.write((char*)C,sigma*sizeof(T));
        fout.write((char*)Occ1[0],sigma*Occ1_sz*sizeof(T));
        fout.write((char*)Occ2[0],sigma*Occ2_sz*sizeof(uint16_t));
        fout.write((char*)SSA,SSA_sz*sizeof(T));
        fout.write((char*)RSRA[0],SSA_sz*sizeof(uint16_t));
        fout.close();
    }
    
    void reindex(BinWord<T> & G, std::string infile)
    {
        std::ifstream fin(infile, std::ios::binary);
        T BV_sz;
        fin.read((char*)&BV_sz,sizeof(T));
        T sz=(BV_sz+7)/8;
        char* bytes=new char[sz];
        fin.read(bytes,sz);
        G.from_chars(bytes,BV_sz);
        fin.read(bytes,sz);
        this->from_chars(bytes,BV_sz);
        delete[] bytes;
        T RSRA_sz, SSA_sz, Occ1_sz, Occ2_sz;
        fin.read((char*)&RSRA_sz,sizeof(T));
        fin.read((char*)&SSA_sz,sizeof(T));
        fin.read((char*)&Occ1_sz,sizeof(T));
        fin.read((char*)&Occ2_sz,sizeof(T));
        Initial((sigma*(Occ1_sz+1)+SSA_sz+3)*sizeof(T)+(sigma+1)*sizeof(T*)+(sigma*Occ2_sz+SSA_sz+2)*sizeof(uint16_t)+(sigma+RSRA_sz+2)*sizeof(uint16_t*));
        C=heap_alloc<T>(sigma);
        Occ1=heap_alloc<T*>(sigma);
        Occ2=heap_alloc<uint16_t*>(sigma);
        RSRA=heap_alloc<uint16_t*>(RSRA_sz);
        RSRA[0]=heap_alloc<uint16_t>(SSA_sz);
        SSA=heap_alloc<T>(SSA_sz);
        Occ1[0]=heap_alloc<T>(sigma*Occ1_sz);
        Occ2[0]=heap_alloc<uint16_t>(sigma*Occ2_sz);
        for(T i=1; i<sigma; ++i)
        {
            Occ1[i]=&Occ1[i-1][Occ1_sz];
            Occ2[i]=&Occ2[i-1][Occ2_sz];
        }
        for(T i=1; i<RSRA_sz; ++i)
        {
            T tmp;
            fin.read((char*)&tmp,sizeof(T));
            RSRA[i]=RSRA[0]+tmp;
        }
        fin.read((char*)C,sigma*sizeof(T));
        fin.read((char*)Occ1[0],sigma*Occ1_sz*sizeof(T));
        fin.read((char*)Occ2[0],sigma*Occ2_sz*sizeof(uint16_t));
        fin.read((char*)SSA,SSA_sz*sizeof(T));
        fin.read((char*)RSRA[0],SSA_sz*sizeof(uint16_t));
        fin.close();
    }
    
    void RelativeOrder(BinWord<T> & G, T* SAI, T* SAIn, T n, T ll)
    {
        T tmp=G((n+1)*ll-1);
        SAIn[2*ll-1]=C[tmp-1]+GetOcc(SAI[0], tmp);
        for(T i=ll-2; i>=0; --i)
        {
            tmp=G(n*ll+i);
            SAIn[i+ll]=C[tmp-1]+GetOcc(SAIn[i+ll+1], tmp);
        }
    }

    std::tuple<T*,T*,T> SegmentSort(BinWord<T> & G, T n, T ll, T* SAmem, T* SAImem)
    {
        T* SA=SAmem+(n%para_sz)*(2*ll-1);
        T* SAI=SAImem+(n%para_sz)*2*ll;
        GetSA(G, n*ll, (n+2)*ll-1, SA, SAI);
        T p=GetALS(G, n*ll, (n+2)*ll-1, SA, SAI, ll);
        return std::make_tuple(SA,SAI,p);
    }

    void CountBWT(T BWT_sz)
    {
        T j1=0;
        for(T i=0, j2=0; i<BWT_sz; ++i)
        {
            if(i%oo2==0)
            {
                ++j1;
                ++j2;
                for(T j=0; j<sigma; ++j)
                {
                    Occ1[j][j1]=Occ1[j][j1-1];
                    Occ2[j][j2]=0;
                }
            }
            else
            {
                if(i%oo==0)
                {
                    ++j2;
                    for(T j=0; j<sigma; ++j)
                        Occ2[j][j2]=Occ2[j][j2-1];
                }
            }
            T tmp=(*this)(i)-1;
            if(tmp>=0)
            {
                ++Occ1[tmp][j1];
                ++Occ2[tmp][j2];
            }
        }
        for(T i=1; i<sigma; ++i)
            C[i]=C[i-1]+Occ1[i-1][j1];
    }
    
    T GetOcc(T i, T c)
    {
        T tmp=Occ1[c-1][(i+1)/oo2];
        if(((i+1)/oo)%oo!=0)
            tmp+=Occ2[c-1][(i+1)/oo];
        for(T j=i-(i+1)%oo+1; j<=i; ++j)
            if((*this)(j)==c)
                ++tmp;
        return tmp;
    }
    
    std::tuple<T,T*> OffSet(T* S, T S_sz)
    {
        T tmp=S[S_sz-1];
        T i=C[tmp-1]+1, j;
        if(tmp==sigma)
            j=this->size()-1;
        else
            j=C[tmp];
        for(T p=S_sz-2; p>=0; --p)
            PreRange(i,j,S[p]);
        T O_sz=j-i+1;
        T* O=new T[O_sz];
        for(T k=i; k<=j; ++k)
            O[k-i]=SimSuffix(k);
        return std::make_tuple(O_sz,O);
    }
    
    T SimSuffix(T i)
    {
        T j=0, k;
        while((k=CheckMember(i))<0)
        {
            T tmp=(*this)(i);
            i=C[tmp-1]+GetOcc(i,tmp);
            ++j;
        }
        return k+j;
    }
    
    T CheckMember(T i)
    {
        uint16_t* p=RSRA[i/oo2];
        uint16_t* q=RSRA[i/oo2+1]-1;
        i=i%oo2;
        if(p>q || *p>i || *q<i)
            return -1;
        else
            if(*p==i)
                return SSA[p-RSRA[0]];
            else
                if(*q==i)
                    return SSA[q-RSRA[0]];
        while(q-p>1)
        {
            uint16_t* m=p+(q-p)/2;
            if(*m==i)
                return SSA[m-RSRA[0]];
            else
                if(*m<i)
                    p=m;
                else
                    q=m;
        }
        return -1;
    }
    
    void PreRange(T& i, T& j, T c)
    {
        i=C[c-1]+GetOcc(i-1,c)+1;
        j=C[c-1]+GetOcc(j,c);
    }
    
    T GetALS(BinWord<T> & G, T h, T t, T* SA, T* SAI, T ll)
    {
        T sz=t-h;
        for(T i=0; i<ll; ++i)
            SA[SAI[i]]=i;
        for(T i=0, j=0; i<sz; ++i)
            if(SA[i]>=0 && SA[i]<ll)
            {
                SA[j]=SA[i];
                SAI[SA[j]]=j;
                ++j;
            }
        for(T i=0, k=0; i<ll; ++i)
            if(SAI[i]!=ll-1)
            {
                T j=SA[SAI[i]+1];
                while(i+k<sz && j+k<sz && G(i+k+h)==G(j+k+h))
                    ++k;
                SAI[SAI[i]+ll]=k;
                if(k>0)
                    --k;
            }
        T p=0;
        for(T i=0; i<ll-1; ++i)
            if(SAI[i+ll]>=ll)
            {
                if(i==0 || SAI[i-1+ll]<ll)
                    SAI[p++]=i;
                if(i==ll-2 || SAI[i+1+ll]<ll)
                    SAI[p++]=i+1;
            }
        return p;
    }

    void GetSA(BinWord<T> & G, T h, T t, T* SA, T* SAI)
    {
        T sz=t-h;
        G.to_T(h,t,SAI);
        for(T i=0; i<sz; ++i)
        {
            SAI[i]+=sz-1;
            SA[i]=i;
        }
        SAI[sz]=-1;
        SortSplit(SA, 0, sz-1, SAI, SAI, true);
        T* SAIh=SAI+this->r;
        while(SA[0]>-sz)
        {
            T i=0;
            T k=0;
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
                    T j=SAI[SA[i]];
                    SortSplit(SA, i, j, SAI, SAIh, true);
                    i=j+1;
                }
            }
            while(i<sz);
            if(k<0)
                SA[i+k]=k;
            SAIh+=(SAIh-SAI);
        }
    }

    void SortSplit(T* SA, T p1, T p2, T* SAI, T* SAIh, bool group)
    {
        T n=p2-p1+1;
        if(n<7)
        {
            InsertionSort(SA, p1, p2, SAI, SAIh, group);
            return;
        }
        T m=p1+n/2;
        T v=SAIh[SA[m]];
        if(n>7)
        {
            T v1=SAIh[SA[p1]], v2=SAIh[SA[p2]];
            if(n>40)
            {
                T s=n/8;
                v1=med3(v1,SAIh[SA[p1+s]],SAIh[SA[p1+2*s]]);
                v=med3(SAIh[SA[m-s]],v,SAIh[SA[m+s]]);
                v2=med3(SAIh[SA[p2-2*s]],SAIh[SA[p2-s]],v2);
            }
            v=med3(v1,v,v2);
        }
        T pa=p1, pb=p1, pc=p2, pd=p2;
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
        T s=std::min(pa-p1,pb-pa);
        for(T i=p1,j=pb-s; i<p1+s; ++i,++j)
            std::swap(SA[i],SA[j]);
        s=std::min(p2-pd,pd-pc);
        for(T i=p2,j=pc+s; i>p2-s; --i,--j)
            std::swap(SA[i],SA[j]);
        pa=p1+pb-pa;
        pc=p2-pd+pc;
        if(pa>p1)
            SortSplit(SA, p1, pa-1, SAI, SAIh, group);
        if(group)
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
            SortSplit(SA, pc+1, p2, SAI, SAIh, group);
    }

    T med3(T a, T b, T c)
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

    void InsertionSort(T* SA, T p1, T p2, T* SAI, T* SAIh, bool group)
    {
        for(T i=p1+1; i<=p2; ++i)
        {
            T j=i;
            while(j>p1 && SAIh[SA[j]]<SAIh[SA[j-1]])
            {
                std::swap(SA[j],SA[j-1]);
                --j;
            }
        }
        if(group)
        {
            T k=p1;
            T v=SAIh[SA[k]];
            for(T i=p1+1; i<=p2; ++i)
            {
                T u=SAIh[SA[i]];
                if(u>v)
                {
                    for(T j=k; j<i; ++j)
                        SAI[SA[j]]=i-1;
                    if(i-k==1)
                        SA[k]=-1;
                    k=i;
                    v=u;
                }
            }
            for(T j=k; j<=p2; ++j)
                SAI[SA[j]]=p2;
            if(k==p2)
                SA[k]=-1;
        }
    }
};

template <typename T>
void index_global(std::string file, T ll, thread_pool & threads1, thread_pool & thread2)
{
    int index_n;
    std::ifstream fin(file);
    fin >> index_n;
    for(int i=0; i<index_n; ++i)
    {
        int file_n;
        fin >> file_n;
        std::list<std::string> files;
        for(int i=0; i<file_n; ++i)
        {
            files.emplace_back();
            fin >> files.back();
        }
        std::string outfile;
        fin >> outfile;
        BinWord<T> G;
        G.readin(files, true);
        BroWheel<T> Bro;
        Bro.index(G, ll, threads1, thread2, outfile);
    }
}

#endif
