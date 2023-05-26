#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "headers/RandomReads.h"

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>

std::condition_variable cv;
std::mutex mtx;
size_t block_num;

void process_and_write_single(std::queue<std::pair<std::string, std::string>> reads, Align &align)
{
    while (!reads.empty())
    {
        align.Oname.swap(reads.front().first);
        align.O.swap(reads.front().second);
        align.Mix();
        align.GetMinimalGraph();
        reads.pop();
    }
    --block_num;
}

void process_and_write(std::queue<std::pair<std::string, std::string>> reads, std::map<std::thread::id, Align> &aligns)
{
    auto &align=aligns.at(std::this_thread::get_id());
    while (!reads.empty())
    {
        align.Oname.swap(reads.front().first);
        align.O.swap(reads.front().second);
        align.Mix();
        align.GetMinimalGraph();
        reads.pop();
    }

    std::unique_lock lck(mtx); // mt
    --block_num;
    cv.notify_all(); // mt
}

std::list<std::future<void>> load_blocks(std::string read_file, thread_pool &threads1, std::map<std::thread::id, Align> &aligns, size_t max_block, size_t block_size)
{
    std::ifstream fin(read_file);
    std::list<std::future<void>> futures;
    do
    {
        std::string name;
        std::string read;
        while (fin >> name >> read)
        {
            std::queue<std::pair<std::string, std::string>> reads;
            reads.emplace(name, read);
            while (reads.size() < block_size && fin >> name >> read)
                reads.emplace(name, read);
            std::unique_lock lck(mtx);
            cv.wait(lck, [max_block] { return block_num < max_block; });
            futures.push_back(threads1.submit(std::bind(process_and_write, std::move(reads), std::ref(aligns))));
            ++block_num;
        }
    } while (fin.good());
    return futures;
}

void load_blocks(std::string read_file, Align &align, size_t block_size)
{
    std::ifstream fin(read_file);
    do
    {
        std::string name;
        std::string read;
        while (fin >> name >> read)
        {
            std::queue<std::pair<std::string, std::string>> reads;
            reads.emplace(name, read);
            while (reads.size() < block_size && fin >> name >> read)
                reads.emplace(name, read);
            ++block_num;
            process_and_write_single(std::move(reads), align);
        }
    } while (fin.good());
}

void test_RandomReads();

int main(int argc, char **argv)
{
    test_RandomReads();
    return 0;
}

void test_RandomReads()
{
    int n_sz = 5, r_sz = 2, t_sz = 1, max_e_sz = 8, lseqlb = 100, lsequb = 200, gseqlb = 10000, gsequb = 20000, aseqlb = 50, asequb = 100, seq_num = 10000, head_in = 10, tail_in = 10;
    std::string argfile("argfile"), local_file("local_file"), global_file("global_file"), read_file("read_file"), truth_file("truth_file");
    bool acyclic = false;
    double gpro = 0.5, rpro = 0, nve=0, nue=0, ve = -5, ue = -2, vf = -5, uf = -2, T = -10, vfp = 0, ufp = 0, vfm = 0, ufm = 0, mat = 1, mis = -3, indel_rate = 0.005, mut_rate = 0.005, apro = 0.5;

    if (false)
    {
        random_DG(n_sz, r_sz, t_sz, max_e_sz, argfile, acyclic, gpro, rpro, nve, nue, ve, ue, vf, uf, T, vfp, ufp, vfm, ufm, mat, mis, local_file, global_file, lseqlb, lsequb, gseqlb, gsequb, aseqlb, asequb, seq_num, read_file, truth_file, indel_rate, mut_rate, head_in, tail_in, apro);
    }

    std::vector<std::string> args;
    std::string tmp;
    std::ifstream fin(argfile);
    while (fin >> tmp)
        args.push_back(tmp.substr(1, tmp.size() - 3));
    fin.close();
    int argc_r = args.size() + 1;
    char **argv_r = new char *[argc_r];
    for (int i = 1; i < argc_r; ++i)
        argv_r[i] = (char *)args[i - 1].c_str();

    std::map<std::string, std::string> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    const std::experimental::filesystem::path current_path(".");
    std::experimental::filesystem::directory_iterator end_itr; // default construction yields past-the-end
    for (std::experimental::filesystem::directory_iterator itr(current_path); itr != end_itr; ++itr)
        if (!std::experimental::filesystem::is_directory(itr->status()))
        {
            std::string file = itr->path().filename().string();
            if (local_file == file.substr(0, local_file.size()))
            {
                std::ifstream fin(file);
                std::string tmp;
                while (fin >> tmp)
                    file2seq[file] += tmp;
            }
            else if (global_file == file.substr(0, global_file.size()))
                file2browheel[file].readin(std::list<std::string>{file}, true);
        }

    int ll = 300, threads1_sz = 10;
    thread_pool threads1(threads1_sz); // mt
    {
        // thread_pool threads1(threads1_sz); // st
        thread_pool thread2(1);
        for (auto &pair : file2browheel)
            pair.second.index(ll, threads1, thread2);
    }

    size_t chunk_sz = size_t(128) * 1024 * 1024;
    std::map<std::thread::id, Align> aligns;
    auto thread_ids=threads1.get_ids();
    std::vector<std::string> mg_files;
    for (int i = 0; i < threads1_sz; ++i)
    {
        mg_files.push_back(read_file + std::to_string(i) + ".mg");
        // aligns.emplace(std::piecewise_construct,std::forward_as_tuple(thread_ids[i]),std::forward_as_tuple(argc_r, argv_r, file2seq, file2browheel, chunk_sz, mg_files.back()));
    }
    // aligns.begin()->second.draw("graph.gv");

    // size_t block_size = 100, max_block = 2 * threads1_sz;
    // std::list<std::future<void>> futures = load_blocks(read_file, threads1, aligns, max_block, block_size); // mt
    // for (auto &future : futures) // mt
    //     future.wait(); // mt
    // // load_blocks(read_file, aligns.begin()->second, block_size); // st

    // for (auto &pair : aligns)
    // {
    //     auto &align=pair.second;
    //     if (!align.Oname.empty())
    //         align.fout.write((char *)&align.max_id, sizeof(align.max_id));
    //     align.fout.close();
    // }

    size_t max_seq = INT64_MAX;
    size_t max_track = 1;
    Track track(argc_r, argv_r, file2seq, file2browheel, chunk_sz);
    track.ReadTrack(mg_files, read_file, read_file + ".trk", read_file + ".alg", max_seq, max_track);
    size_t max_extract = max_track;
    track.extract(read_file + ".trk", read_file + ".ext", max_extract);
}