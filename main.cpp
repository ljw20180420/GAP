#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>
#include "headers/test.h"

typedef ssize_t temp_type;
size_t block_num;
std::mutex mtx;
std::condition_variable cv;

void process_and_write(std::queue<std::string> reads, std::map<std::thread::id,Align<temp_type>> & aligns, size_t block_index, std::string min_graph_file)
{
    Align<temp_type> & align=aligns.at(std::this_thread::get_id());
//     Align<temp_type> & align=aligns.begin()->second; // debug
    align.MGout.open(min_graph_file+std::to_string(block_index));
    while(!reads.empty())
    {
        for(size_t i=0; i<reads.front().size(); ++i)
            reads.front()[i]=(reads.front()[i]/2+3) & 7;
        align.Mix(reads.front());
        align.CallGetMinimalGraph(reads.front().size(),reads.front());
        reads.pop();
    }
    align.MGout.close();
    
    std::unique_lock lck(mtx); // mt
    cv.notify_all(); // mt
    --block_num;
}

bool load_blocks(std::ifstream & fin, std::queue<std::future<void>> & futures, thread_pool & threads1, std::map<std::thread::id,Align<temp_type>> & aligns, size_t max_block, size_t block_size, size_t & block_index, std::string & min_graph_file)
{
    std::string read;
    while(block_num<max_block && fin>>read)
    {
        std::queue<std::string> reads;
        reads.push(read);
        while(reads.size()<block_size && fin>>read)
            reads.push(read);
        futures.push(threads1.submit(std::bind(process_and_write, std::move(reads), std::ref(aligns), block_index, min_graph_file)));
//         process_and_write(std::move(reads), std::ref(aligns), block_index, min_graph_file); // debug
        ++block_num;
        ++block_index;
    }
    fin>>std::ws;
    return fin.good();
}

void load_and_distribute(std::string read_file, std::string min_graph_file, std::string path_file, std::string predict_file, size_t max_block, size_t block_size, thread_pool & threads1, std::string graph_file, size_t chunk_sz, size_t chunk_sz_2, int max_num, int max_glo)
{
    to_gv(graph_file);
    Graph<temp_type> graph(graph_file);
    std::queue<std::future<void>> futures;
    std::map<std::thread::id,Align<temp_type>> aligns;
    std::vector<std::thread::id> thread_ids=threads1.get_ids();
    for(size_t i=0; i<thread_ids.size(); ++i)
        aligns.emplace(std::piecewise_construct,std::forward_as_tuple(thread_ids[i]),std::forward_as_tuple(graph,chunk_sz,chunk_sz_2));
    std::ifstream fin(read_file);
    block_num=0;
    size_t block_index=0;
    std::unique_lock lck(mtx); // mt
    while(load_blocks(fin,futures,threads1,aligns,max_block,block_size,block_index,min_graph_file))
        cv.wait(lck); // mt
//         ;
    lck.unlock(); // mt
    fin.close();
    
    std::ofstream fout(path_file);
    fout << std::setprecision(10);
    size_t read_num=0;
    Tracker<temp_type> tracker;
    for(size_t bi=0; bi<block_index; ++bi)
    {
        futures.front().wait();
        futures.pop();
        tracker.ReadAndExploreAll(min_graph_file,bi,fout,max_num,max_glo,graph,read_num);
    }
    fout.close();
    Truther<temp_type> truther;
    truther.PathToTruth(path_file,predict_file,max_glo);
}

int main(int argc, char **argv) 
{
    // warning if node with nonzero ue has global out edge
    bool acyclic=false;
    int n_sz=5, r_sz=2, lseqlb=100, lsequb=200, gseqlb=10000, gsequb=20000, aseqlb=50, asequb=100, seq_num=10000, head_in=10, tail_in=10;
    double gpro=0.5, rpro=0, tve=-(1e-6), tue=0, ve=-5, ue=-2, vf=-5, uf=-2, T=-10, tvf=-(1e-6), tuf=0, mat=1, mis=-3, indel_rate=0.005, mut_rate=0.005, apro=0.5;
    
    int ll=300, max_num=1, max_glo=3, threads1_sz=10;
    size_t kb=1024;
    size_t chunk_sz=size_t(128)*kb*kb, chunk_sz_2=size_t(128)*kb*kb, block_size=100, max_block=2*threads1_sz;
    std::string graph_file("graph_file"), local_file("local_file"), local_file_score("local_file_score"), index_file("index_file"), global_file("global_file"), global_file_score("global_file_score"), index_information("index_information"), read_file("read_file"), truth_file("truth_file"), min_graph_file("tmp/min_graph_file"), path_file("path_file"), predict_file("predict_file");
    
    random_DG(n_sz, r_sz, gpro, rpro, tve, tue, graph_file, local_file, local_file_score, index_file, global_file, global_file_score, lseqlb, lsequb, gseqlb, gsequb, index_information, ve, ue, vf, uf, T, tvf, tuf, mat, mis, aseqlb, asequb, seq_num, read_file, truth_file, indel_rate, mut_rate, head_in, tail_in, acyclic, apro);
    
    thread_pool threads1(threads1_sz);
    {
        thread_pool thread2(1);
        index_global<temp_type>(index_information,ll,threads1,thread2);
    }
    
    load_and_distribute(read_file,min_graph_file,path_file,predict_file,max_block,block_size,threads1,graph_file,chunk_sz,chunk_sz_2,max_num,max_glo);
    
    return 0;
}
