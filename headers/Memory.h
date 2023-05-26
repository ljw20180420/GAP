#ifndef MEMORY_H
#define MEMORY_H

#include <vector>

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

#endif