#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "headers/RandomReads.h"

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>
#include <unistd.h>
#include <ctime>

struct Parallel_Align
{
    int argc;
    char **argv;
    std::string run_name;
    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    int Omax;
    std::deque<std::string> read_files, mg_files;
    std::atomic<int> block_num;

    Parallel_Align(int argc_, char **argv_, std::string run_name_)
        : argc(argc_), argv(argv_), run_name(run_name_)
    {
        initialize();
    }

    Parallel_Align(std::string argfile, std::string run_name_)
        : run_name(run_name_)
    {
        std::vector<std::string> args;
        std::string tmp;
        std::ifstream fin(argfile);
        while (fin >> tmp)
        {
            if (tmp.back() == ',')
                args.push_back(tmp.substr(1, tmp.size() - 3));
            else
                args.push_back(tmp.substr(1, tmp.size() - 2));
        }
        fin.close();
        argc = args.size() + 1;
        argv = new char *[argc];
        for (int i = 1; i < argc; ++i)
        {
            argv[i] = new char[args[i - 1].size() + 1];
            strcpy(argv[i], args[i - 1].c_str());
        }
        initialize();
    }

    void initialize()
    {
        Graph(argc, argv, file2seq, file2browheel).draw("graph.gv");
        for (int i = 1; i < argc; ++i)
            if (!strcmp(argv[i], "--read_files"))
            {
                while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                    read_files.push_back(argv[i]);
                --i;
            }
        Omax = 0;
        for (auto &read_file : read_files)
        {
            std::ifstream fin(read_file);
            int maxlen = 0;
            int64_t preg = 0;
            while (fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
            {
                int64_t nowg = fin.tellg();
                int len = nowg - preg;
                if (len > maxlen)
                    maxlen = len;
                preg = nowg;
            }
            --maxlen;
            if (maxlen > Omax)
                Omax = maxlen;
        }
    }

    void do_index(int ll, int threads1_sz, bool reverse_complement)
    {
        for (auto &pair : file2browheel)
            pair.second.readin(pair.first, true, reverse_complement);
        thread_pool threads1(threads1_sz), thread2(1);
        for (auto &pair : file2browheel)
        {
            pair.second.index(ll, threads1, thread2);
            pair.second.saveBroWheel();
        }
    }

    void load_index()
    {
        for (auto &pair : file2seq)
            pair.second.readin(pair.first);
        for (auto &pair : file2browheel)
            pair.second.loadBroWheel(pair.first);
    }

    std::queue<std::pair<std::string, std::string>> fill_block(int64_t block_size, std::ifstream &fin, std::deque<std::string>::iterator &iter)
    {
        std::queue<std::pair<std::string, std::string>> reads;
        while (reads.size() < block_size)
        {
            std::string name, read;
            if (std::getline(std::getline(fin, name), read))
            {
                reads.emplace(name, read);
                if (name[0] == '@')
                    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n').ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            else
            {
                fin.close();
                if (++iter != read_files.end())
                    fin.open(*iter);
                else
                    break;
            }
        }
        return reads;
    }

    void process_and_write_single(std::queue<std::pair<std::string, std::string>> &&reads, Align &align)
    {
        while (!reads.empty())
        {
            align.Oname.swap(reads.front().first);
            align.O.swap(reads.front().second);
            align.Mix();
            align.GetMinimalGraph();
            reads.pop();
        }
    }

    void single_align(int64_t block_size)
    {
        mg_files.emplace_back(run_name + ".mg");
        Align align(argc, argv, file2seq, file2browheel, mg_files.front(), Omax);

        auto iter = read_files.begin();
        std::ifstream fin(*iter);
        while (iter != read_files.end())
        {
            process_and_write_single(fill_block(block_size, fin, iter), align);
        }
        if (!align.Oname.empty())
            align.fout.write((char *)&align.max_id, sizeof(align.max_id));
        align.fout.close();
    }

    void parallel_align(int threads1_sz, int64_t block_size)
    {
        std::map<std::thread::id, Align> aligns;
        thread_pool threads1(threads1_sz);
        auto thread_ids = threads1.get_ids();
        for (size_t i = 0; i < threads1.size(); ++i)
        {
            mg_files.emplace_back(run_name + std::to_string(i) + ".mg");
            aligns.emplace(std::piecewise_construct, std::forward_as_tuple(thread_ids[i]), std::forward_as_tuple(argc, argv, file2seq, file2browheel, mg_files.back(), Omax));
        }

        int max_block = 2 * threads1_sz;
        block_num.store(0);
        std::deque<std::future<void>> futures;
        auto iter = read_files.begin();
        std::ifstream fin(*iter);
        while (iter != read_files.end())
        {
            auto reads = fill_block(block_size, fin, iter);
            while (block_num.load() == max_block)
                ;
            futures.push_back(threads1.submit(std::bind(&Parallel_Align::process_and_write, this, std::move(reads), std::ref(aligns))));
            ++block_num;
        }
        for (auto &future : futures)
            future.wait();

        for (auto &pair : aligns)
        {
            auto &align = pair.second;
            if (!align.Oname.empty())
                align.fout.write((char *)&align.max_id, sizeof(align.max_id));
            align.fout.close();
        }
    }

    void process_and_write(std::queue<std::pair<std::string, std::string>> reads, std::map<std::thread::id, Align> &aligns)
    {
        process_and_write_single(std::move(reads), aligns.at(std::this_thread::get_id()));
        --block_num;
    }

    void track(int64_t max_seq, int64_t max_track, int64_t max_extract)
    {
        Track track(argc, argv, file2seq, file2browheel);
        track.ReadTrack(mg_files, read_files, run_name, max_seq, max_track, max_extract, 0);
    }
};

void help()
{
    std::cout << "--read_files file1 file2 ...\n"
              << "--nodes node_name1 gap_open_penalty(GOP)1 gap_extend_penalty(GEP)1 node_name2...\n"
              << "--roots root_name1 root_name2 ...\n"
              << "--targets target_name1 target_name2 ...\n"
              << "--locals <default_gamma | 49_pair_scores1> local_file1 GOP_in_ref1 GEP_in_ref1 GOP_in_query1 GEP_in_query1 basic_penalty1 GOP_in_ref_begin1 GEP_in_ref_begin1 GOP_in_ref_end1 GEP_in_ref_end1 tail1 head1 <default_gamma | 49_pair_scores2> local_file2...\n"
              << "--globals <default_gamma | 49_pair_scores1> global_file1 GOP_in_ref1 GEP_in_ref1 GOP_in_query1 GEP_in_query1 basic_penalty1 tail1 head1 <default_gamma | 49_pair_scores2> global_file2...\n";
}

std::deque<std::string> get_dirs()
{
    const std::experimental::filesystem::path current_path(".");
    std::experimental::filesystem::directory_iterator end_itr; // default construction yields past-the-end
    std::deque<std::string> dirs;
    for (std::experimental::filesystem::directory_iterator itr(current_path); itr != end_itr; ++itr)
        if (std::experimental::filesystem::is_directory(itr->status()))
            dirs.push_back(itr->path().filename().string());
    return dirs;
}

std::deque<std::string> get_argfiles()
{
    const std::experimental::filesystem::path current_path(".");
    std::experimental::filesystem::directory_iterator end_itr; // default construction yields past-the-end
    std::deque<std::string> argfiles;
    for (std::experimental::filesystem::directory_iterator itr(current_path); itr != end_itr; ++itr)
        if (!std::experimental::filesystem::is_directory(itr->status()))
        {
            std::string str_tmp = itr->path().filename().string();
            if (str_tmp.size() > 8 && str_tmp.substr(str_tmp.size() - 8, 8) == ".argfile")
                argfiles.push_back(str_tmp);
        }
    return argfiles;
}

void test_GenerateRandomReads(std::string dir, std::string argfile)
{
    chdir(dir.c_str());

    int n_sz = 5, r_sz = 2, t_sz = 1, max_e_sz = 8, lseqlb = 100, lsequb = 200, gseqlb = 10000, gsequb = 20000, aseqlb = 50, asequb = 100, seq_num = 10000, head_in = 10, tail_in = 10;
    bool acyclic = false;
    double gpro = 0.5, rpro = 0, nve = 0, nue = 0, ve = -5, ue = -2, vf = -5, uf = -2, T = -10, dT = -5, minScore = 20, vfp = 0, ufp = 0, vfm = 0, ufm = 0, mat = 1, mis = -3, indel_rate = 0.005, mut_rate = 0.005, apro = 0.5;
    int64_t diffseg = 10;

    random_DG(n_sz, r_sz, t_sz, max_e_sz, "argfile", acyclic, gpro, rpro, nve, nue, ve, ue, vf, uf, T, dT, minScore, diffseg, vfp, ufp, vfm, ufm, mat, mis, "local_file", "global_file", lseqlb, lsequb, gseqlb, gsequb, aseqlb, asequb, seq_num, "read_file", "truth_file", indel_rate, mut_rate, head_in, tail_in, apro);
}

void test_Align(std::string dir, std::string run_name, std::string argfile, int threads1_sz, int64_t max_track, int64_t max_extract)
{
    chdir(dir.c_str());

    Parallel_Align parallel_align(argfile, run_name);
    // parallel_align.do_index(150000000, 3, true);  
    parallel_align.load_index();
    // time_t time_tic = time(0);
    // // parallel_align.single_align(100); // st
    // parallel_align.parallel_align(threads1_sz, 100); // mt
    // std::cerr << "align time is " << time(0) - time_tic << '\n';

    for (size_t i = 0; i < threads1_sz; ++i) // test track
        parallel_align.mg_files.emplace_back(run_name + std::to_string(i) + ".mg"); // test track

    parallel_align.track(INT64_MAX, max_track, max_extract);
}

void run_sim(std::string sim_dir, std::string run_name, int threads1_sz, int64_t max_track, int64_t max_extract)
{
    chdir(sim_dir.c_str());
    std::deque<std::string> dirs = get_dirs();
    for (auto dir : dirs)
    {
        chdir(("./" + dir).c_str());
        std::deque<std::string> argfiles = get_argfiles();
        for (auto &argfile : argfiles)
            test_Align("./", run_name, argfile, threads1_sz, max_track, max_extract);
        chdir("..");
    }
    chdir("..");
}

int main(int argc, char **argv)
{
    // test_GenerateRandomReads("./test_RandomReads", "argfile");
    test_Align("./test_RandomReads", "random", "argfile", 24, 1, 1);

    
    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time/test_CRISPR_5", "CRISPR", "argfile", 6, 1, 1);
    
    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time/test_CRISPR_10", "CRISPR", "argfile", 6, 1, 1);

    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time/test_CRISPR_15", "CRISPR", "argfile", 6, 1, 1);

    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time/test_CRISPR_20", "CRISPR", "argfile", 6, 1, 1);

    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time/test_CRISPR_25", "CRISPR", "argfile", 6, 1, 1);

    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_CRISPR_cross", "CRISPR", "argfile", 6, 1, 1);

    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/zyj/zhengjin/4C_fold/run", "4C", "argfile", 12, 1, 1);

    // run_sim("./test_single_cut", "sim", 12, 10, 10);
    // run_sim("./test_double_cut", "sim", 12, 10, 10);

    return 0;
}