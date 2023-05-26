#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "headers/BroWheel.h"

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>

int main(int argc, char **argv) 
{
    BroWheel<ssize_t> browheel;
    std::list<std::string>  files{"global_file0","global_file1","global_file2"};
    browheel.readin(files, false);
    ssize_t ll=100;
    thread_pool threads1(6);
    thread_pool thread2(1);
    browheel.index(ll, threads1, thread2);
    browheel.saveBroWheel("index");
    BroWheel<ssize_t> browheel2;
    browheel2.loadBroWheel("index");
    std::vector<ssize_t> seq(browheel2.sequence.size());
    std::ofstream fout("sequence");
    for(ssize_t i=0; i<browheel2.sequence.size(); ++i)
    {
        seq[i]=browheel2.sequence(i);
        fout << seq[i];
    }
    fout.close();

    for(ssize_t s=0; s<browheel2.sequence.size(); s+=999)
        for(ssize_t l=1000; l<=browheel2.sequence.size()-s; l+=1000)
        {
            ssize_t occ_sz;
            ssize_t* occ;
            std::tie(occ_sz,occ)=browheel2.OffSet(seq.data()+s, l);
            for(ssize_t o=0; o<occ_sz; ++o)
            {
                for(ssize_t i=0; i<l; ++i)
                    if(seq[occ[o]+i]!=seq[s+i])
                        std::cerr << "error" << std::endl;
            }
        }
    
    return 0;
}
