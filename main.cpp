#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>
// #include "headers/Threadpool.h"
// #include "headers/BroWheel.h"
// #include "headers/Graph.h"
#include "headers/test.h"

const int threads1_sz=6;
typedef int temp_type;

int main(int argc, char **argv) 
{
    // warning if node with nonzero ue has global out edge
    bool acyclic=false;
    int n_sz=5, r_sz=2, lseqlb=100, lsequb=200, gseqlb=10000, gsequb=20000, aseqlb=50, asequb=100, seq_num=100, head_in=10, tail_in=10, ll=300, max_num=1, max_glo=3;
    double gpro=0.5, rpro=0, tve=-(1e-6), tue=0, ve=-5, ue=-2, vf=-5, uf=-2, T=-10, tvf=-(1e-6), tuf=0, mat=1, mis=-3, indel_rate=0.005, mut_rate=0.005, apro=0.5;
    std::string graph_file("graph_file"), local_file("local_file"), local_file_score("local_file_score"), index_file("index_file"), global_file("global_file"), global_file_score("global_file_score"), index_information("index_information"), read_file("read_file"), truth_file("truth_file"), min_graph_file("min_graph_file"), path_file("path_file"), predict_file("predict_file");
    random_DG(n_sz, r_sz, gpro, rpro, tve, tue, graph_file, local_file, local_file_score, index_file, global_file, global_file_score, lseqlb, lsequb, gseqlb, gsequb, index_information, ve, ue, vf, uf, T, tvf, tuf, mat, mis, aseqlb, asequb, seq_num, read_file, truth_file, indel_rate, mut_rate, head_in, tail_in, acyclic, apro);
    {
        to_gv(graph_file);
        thread_pool threads1(threads1_sz);
        thread_pool thread2(1);
        index_global<temp_type>(index_information,ll,threads1,thread2);
    }
    Graph<temp_type> graph(graph_file);
    size_t chunk_sz=1024*1024*1024;
    Align<temp_type> align(graph,chunk_sz);
    align.MGout.open(min_graph_file);
    std::ifstream fin(read_file);
    std::string O;
    while(fin >> O)
    {
        for(size_t i=0; i<O.size(); ++i)
            O[i]=(O[i]/2+3) & 7;
        align.Mix(O);
        align.CallGetMinimalGraph(O.size(),O);
    }
    fin.close();
    align.MGout.close();
    Tracker<temp_type> tracker;
    tracker.ReadAndExploreAll(min_graph_file,path_file,max_num,max_glo,graph);
    Truther<temp_type> truther;
    truther.PathToTruth(path_file,predict_file,max_glo);
    return 0;
}
