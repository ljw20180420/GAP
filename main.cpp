#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "headers/BroWheel.h"

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>

const std::string BroWheel::int2base="KLNACTG";

int main(int argc, char **argv) 
{
    BroWheel browheel;
    std::list<std::string>  files{"global_file0","global_file1","global_file2"};
    browheel.readin(files, false);
    size_t ll=100;
    thread_pool threads1(6);
    thread_pool thread2(1);
    browheel.index(ll, threads1, thread2);
    browheel.saveBroWheel("index");
    BroWheel browheel2;
    browheel2.loadBroWheel("index");
    std::vector<int> seq(browheel2.sequence.size());
    std::ofstream fout("sequence");
    for(size_t i=0; i<browheel2.sequence.size(); ++i)
    {
        seq[i]=browheel2.sequence(i);
        fout << seq[i];
    }
    fout.close();

    for(size_t s=0; s<browheel2.sequence.size(); s+=999)
        for(size_t l=1000; l<=browheel2.sequence.size()-s; l+=1000)
        {
            size_t occ_sz;
            size_t* occ;
            std::tie(occ_sz,occ)=browheel2.OffSet(seq.data()+s, l);
            for(size_t o=0; o<occ_sz; ++o)
            {
                for(size_t i=0; i<l; ++i)
                    if(seq[occ[o]+i]!=seq[s+i])
                        std::cerr << "error" << std::endl;
            }
        }
    
    return 0;
}
