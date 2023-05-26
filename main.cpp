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
#include "headers/Help.h"

struct Parallel_Align
{
    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    Graph graph;

    std::string argfile, run_name;
    std::deque<std::string> read_files, mg_files;
    int index_threads_sz, align_threads_sz, max_round, max_extract, diff_thres, max_range, min_seg_num, max_seg_num, block_size, Omax, round = 0;
    int64_t max_mega, ll;
    std::atomic<int> block_num;

    Parallel_Align(std::string argfile_)
        : graph(argfile_, file2seq, file2browheel), argfile(argfile_)
    {
        std::string data;
        std::ifstream fin(argfile);
        std::getline(fin, data, '\0');
        std::map<std::string, std::string> hts = get_handles(split_string(data, ';', '{', '}'));
        run_name = hts["run_name"];
        read_files = split_string(remove_open_close(hts["read_files"], '{', '}'), ';', '\0', '\0');
        index_threads_sz = std::stoi(hts["index_threads_sz"]);
        align_threads_sz = std::stoi(hts["align_threads_sz"]);
        max_round = std::stoi(hts["max_round"]);
        max_extract = std::stoi(hts["max_extract"]);
        diff_thres = std::stoi(hts["diff_thres"]);
        max_range = std::stoi(hts["max_range"]);
        min_seg_num = std::stoi(hts["min_seg_num"]);
        max_seg_num = std::stoi(hts["max_seg_num"]);
        max_mega = std::stol(hts["max_mega"]);
        block_size = std::stoi(hts["block_size"]);
        ll = std::stol(hts["ll"]);
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
        graph.draw("graph.gv");
    }

    void do_index()
    {
        thread_pool threads1(index_threads_sz), thread2(1);
        for (auto &pair : file2browheel)
        {
            int true_num = 0, false_num = 0;
            for (EdgeGlobalCross &global : graph.global_crosses)
                if (global.name == pair.first)
                    if (global.reverse_complement)
                        ++true_num;
                    else
                        ++false_num;
            for (EdgeGlobalCircuit &global : graph.global_circuits)
                if (global.name == pair.first)
                    if (global.reverse_complement)
                        ++true_num;
                    else
                        ++false_num;
            bool reverse_complement;
            if (true_num > 0 && false_num > 0)
            {
                std::cerr << "unable to determine whether " << pair.first << " should be reverse_complement\n";
                exit(-1);
            }
            else if (true_num == 0 && false_num == 0)
                continue;
            else if (true_num > 0)
                reverse_complement = true;
            else
                reverse_complement = false;
            pair.second.readin(pair.first, true, reverse_complement);
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

    std::queue<std::pair<std::string, std::string>> fill_block(std::ifstream &fin, std::deque<std::string>::iterator &iter)
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

    void single_align()
    {
        if (round > 0)
            read_files = {run_name + std::to_string(round - 1) + ".fail"};
        mg_files = {run_name + std::to_string(round) + ".mg"};
        Align align(argfile, file2seq, file2browheel, mg_files.front(), Omax);
        for (Edge *edge : align.edges)
            edge->T += round * edge->dT;

        auto iter = read_files.begin();
        std::ifstream fin(*iter);
        while (iter != read_files.end())
        {
            process_and_write_single(fill_block(fin, iter), align);
        }
        if (!align.Oname.empty())
            align.fout.write((char *)&align.max_id, sizeof(align.max_id));
        align.fout.close();
    }

    void parallel_align()
    {
        if (round > 0)
            read_files = {run_name + std::to_string(round - 1) + ".fail"};
        std::map<std::thread::id, Align> aligns;
        thread_pool threads1(align_threads_sz);
        auto thread_ids = threads1.get_ids();
        mg_files.clear();
        for (size_t i = 0; i < threads1.size(); ++i)
        {
            mg_files.emplace_back(run_name + std::to_string(round) + ".mg" + std::to_string(i));
            aligns.emplace(std::piecewise_construct, std::forward_as_tuple(thread_ids[i]), std::forward_as_tuple(argfile, file2seq, file2browheel, mg_files.back(), Omax));
        }
        for (auto &pair : aligns)
            for (Edge *edge : pair.second.edges)
                edge->T += round * edge->dT;

        int max_block = 2 * align_threads_sz;
        block_num.store(0);
        std::deque<std::future<void>> futures;
        auto iter = read_files.begin();
        std::ifstream fin(*iter);
        while (iter != read_files.end())
        {
            auto reads = fill_block(fin, iter);
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

    int64_t track()
    {
        Track track(argfile, file2seq, file2browheel, run_name + std::to_string(round), max_extract, max_mega, diff_thres, max_range, min_seg_num, max_seg_num);
        return track.ReadTrack(mg_files, read_files);
    }
};

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
    double gpro = 0.5, rpro = 0, nve = 0, nue = 0, ve = -5, ue = -2, vf = -5, uf = -2, T = -10, dT = -5, min_score = 20, vfp = 0, ufp = 0, vfm = 0, ufm = 0, mat = 1, mis = -3, indel_rate = 0.005, mut_rate = 0.005, apro = 0.5;

    int index_threads_sz = 3, align_threads_sz = 24, max_extract = 3, diff_thres = 2, max_range = 10, min_seg_num = 0, max_seg_num = 0, max_round = 5, block_size = 100;
    int64_t max_mega = 10000, ll = 150000000;

    random_DG(n_sz, r_sz, t_sz, max_e_sz, "argfile", acyclic, gpro, rpro, nve, nue, ve, ue, vf, uf, T, dT, min_score, vfp, ufp, vfm, ufm, mat, mis, "local_file", "global_file", lseqlb, lsequb, gseqlb, gsequb, aseqlb, asequb, seq_num, "read_file", "truth_file", indel_rate, mut_rate, head_in, tail_in, apro, "random", index_threads_sz, align_threads_sz, max_extract, diff_thres, max_range, min_seg_num, max_seg_num, max_round, block_size, max_mega, ll);
}

void override(std::experimental::filesystem::path dir, std::string argfile, std::deque<std::pair<std::string, std::string>> hvps)
{
    chdir(dir.c_str());

    std::string data;
    std::ifstream fin(argfile);
    std::getline(fin, data, '\0');
    for (std::pair<std::string, std::string> &hvp : hvps)
    {
        size_t pos = 0;
        while (true)
        {
            pos = data.find(hvp.first, pos);
            if (pos == std::string::npos)
                break;
            else
            {
                int cpos = data.rfind(";", pos);
                if (cpos == std::string::npos)
                    cpos = -1;
                for (++cpos; isspace(data[cpos]); ++cpos)
                    ;
                int epos = data.find("=", pos), eepos = epos;
                if (epos == std::string::npos)
                {
                    pos += hvp.first.size();
                    continue;
                }
                for (--epos; isspace(data[epos]); --epos)
                    ;
                if (cpos < pos || epos > pos + hvp.first.size() - 1)
                {
                    pos += hvp.first.size();
                    continue;
                }
                int depth = 0;
                for (pos = eepos + 1; pos < data.size(); ++pos)
                {
                    if (data[pos] == ';')
                    {
                        if (depth == 0)
                            break;
                    }
                    else if (data[pos] == '{')
                        ++depth;
                    else if (data[pos] == '}')
                        --depth;
                }
                data.replace(eepos + 1, pos - eepos - 1, " " + hvp.second);
                pos += hvp.second.size() + 1 - (pos - eepos - 1);
            }
        }
    }
    std::ofstream fout(argfile + ".override");
    fout << data;
}

std::deque<int64_t> test_Align(std::experimental::filesystem::path dir, std::string argfile, bool reindex, bool para)
{
    chdir(dir.c_str());

    int64_t not_failed_by_unalign;
    Parallel_Align parallel_align(argfile);
    if (reindex)
        parallel_align.do_index();
    else
        parallel_align.load_index();
    std::deque<int64_t> align_times;
    do
    {
        time_t time_tic = time(0);
        if (para)
            parallel_align.parallel_align();
        else
            parallel_align.single_align();
        
        align_times.push_back(time(0) - time_tic);
        std::cerr << "align time is " << align_times.back() << '\n';

        // parallel_align.mg_files.clear(); // test track
        // for (size_t i = 0; i < parallel_align.align_threads_sz; ++i) // test track
        //     parallel_align.mg_files.emplace_back(parallel_align.run_name + std::to_string(parallel_align.round) + ".mg" + std::to_string(i)); // test track
        // if (parallel_align.round > 0) // test track
        //     parallel_align.read_files = {parallel_align.run_name + std::to_string(parallel_align.round - 1) + ".fail"}; // test track

        not_failed_by_unalign = parallel_align.track();

        ++parallel_align.round;

    } while (not_failed_by_unalign > 0 && parallel_align.round <= parallel_align.max_round);

    return align_times;
}

void run_sim(std::string sim_dir, bool para)
{
    chdir(sim_dir.c_str());
    std::deque<std::string> dirs = get_dirs();
    for (auto dir : dirs)
    {
        chdir(("./" + dir).c_str());
        std::deque<std::string> argfiles = get_argfiles();
        for (auto &argfile : argfiles)
            test_Align("./", argfile, false, para);
        chdir("..");
    }
    chdir("..");
}

int gradual_align(std::experimental::filesystem::path dir, std::string argfile, std::string argfile_last, std::string run_name, int align_trheads_sz, int max_extract, int max_mega, int diff_thres, int max_range, int min_seg_num, int max_seg_num, int block_size, int max_round, double T, double dT, int min_score, bool reindex, bool para)
{
    override(dir, argfile, {{"run_name", run_name.c_str()}, {"align_threads_sz", std::to_string(align_trheads_sz).c_str()}, {"max_extract", std::to_string(max_extract).c_str()}, {"max_mega", std::to_string(max_mega).c_str()}, {"diff_thres", std::to_string(diff_thres).c_str()}, {"max_range", std::to_string(max_range).c_str()}, {"min_seg_num", std::to_string(min_seg_num).c_str()}, {"max_seg_num", std::to_string(max_seg_num).c_str()}, {"block_size", std::to_string(block_size).c_str()}, {"max_round", std::to_string(max_round).c_str()}, {"T", std::to_string(T).c_str()}, {"dT", std::to_string(dT).c_str()}, {"min_score", std::to_string(min_score).c_str()}});
    std::deque<int64_t> align_times = test_Align(dir, argfile + ".override", reindex, para);
    std::ofstream fout((dir / "align_times").string());
    for (int round = 0; round < align_times.size(); ++round)
        fout << round << '\t' << T + round * dT << '\t' << align_times[round] << '\n';
    for (int i = 5; i >= 0; --i)
    {
        std::string fail_last = run_name + std::to_string(i) + ".fail";
        std::ifstream fin(fail_last);
        if (!fin)
            continue;
        fin.seekg(0, std::ios::end);
        if (fin.tellg() > 0)
        {
            fin.close();
            override(dir, argfile_last, {{"run_name", (run_name + "_last").c_str()}, {"read_files", "{" + fail_last + ";};"}, {"align_threads_sz", std::to_string(align_trheads_sz).c_str()}, {"max_extract", std::to_string(max_extract).c_str()}, {"max_mega", std::to_string(max_mega).c_str()}, {"diff_thres", std::to_string(diff_thres).c_str()}, {"max_range", std::to_string(max_range).c_str()}, {"min_seg_num", "0"}, {"max_seg_num", "0"}, {"block_size", std::to_string(block_size).c_str()}, {"max_round", "0"}, {"T", "0"}, {"dT", "0"}, {"min_score", "0"}});
            std::deque<int64_t> align_times = test_Align(dir, argfile_last + ".override", false, para);
            fout << "last\t*\t" << align_times.front() << '\n';
        }
        break;
    }
}

int main(int argc, char **argv)
{
    // help();

    // test_GenerateRandomReads("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_RandomReads", "argfile");
    // override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_RandomReads", "argfile", {{"align_threads_sz", "24"}, {"max_extract", "4"}, {"max_mega", "20000"}, {"diff_thres", "1"}, {"max_range", "5"}, {"max_round", "5"}, {"T", "-5"}, {"dT", "-1"}, {"min_score", "20"}});
    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_RandomReads", "argfile.override", false, true);

    // override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_CRISPR", "argfile", {{"align_threads_sz", "6"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"max_round", "4"}, {"T", "-10"}, {"dT", "-5"}, {"min_score", "30"}});
    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_CRISPR", "argfile.override", false, true);

    // override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_CRISPR_cross", "argfile", {{"align_threads_sz", "6"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"max_round", "4"}, {"T", "-10"}, {"dT", "-5"}, {"min_score", "20"}});
    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_CRISPR_cross", "argfile.override", false, true);

    // std::ofstream fout("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time/align_times");
    // for (int i = 0; i < 5; ++i)
    // {
    //     override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time", "argfile", {{"run_name", "CRISPR" + std::to_string(5 * (i + 1))}, {"align_threads_sz", "6"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"min_seg_num", "1"}, {"max_seg_num", "1"}, {"block_size", "10"}, {"max_round", "0"}, {"T", std::to_string(-5 * (i + 1))}, {"min_score", "0"}});
    //     std::deque<int64_t> align_times = test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time", "argfile.override", false, true);
    //     fout << align_times.front() << '\n';
    // }

    // {
    //     std::ofstream fout("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_direct/align_times"), fout2("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_direct/total_times");
    //     override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_direct", "argfile", {{"run_name", "CRISPR"}, {"align_threads_sz", "6"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"min_seg_num", "0"}, {"max_seg_num", "0"}, {"block_size", "10"}, {"max_round", "0"}, {"T", "0"}, {"dT", "0"}, {"min_score", "0"}});
    //     time_t time_tic = time(0);
    //     std::deque<int64_t> align_times = test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_direct", "argfile.override", false, true);
    //     fout << align_times.front() << '\n';
    //     fout2 << time(0) - time_tic << '\n';
    // }

    // {
    //     std::ofstream fout("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop/align_times");
    //     override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop", "argfile", {{"run_name", "CRISPR"}, {"align_threads_sz", "6"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"min_seg_num", "1"}, {"max_seg_num", "1"}, {"block_size", "10"}, {"max_round", "5"}, {"T", "-5"}, {"dT", "-5"}, {"min_score", "0"}});
    //     time_t time_tic = time(0);
    //     std::deque<int64_t> align_times = test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop", "argfile.override", false, true);
    //     for (int64_t align_time : align_times)
    //         fout << align_time << '\n';
    //     for (int i = 5; i >= 0; --i)
    //     {
    //         std::string fail_last = "CRISPR" + std::to_string(i) + ".fail";
    //         std::ifstream fin(fail_last);
    //         if (!fin)
    //             continue;
    //         fin.seekg(0, std::ios::end);
    //         if (fin.tellg() > 0)
    //         {
    //             fin.close();
    //             override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop", "argfile_last", {{"run_name", "CRISPR_last"}, {"read_files", "{" + fail_last + ";};"}, {"align_threads_sz", "6"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"min_seg_num", "0"}, {"max_seg_num", "0"}, {"block_size", "10"}, {"max_round", "0"}, {"T", "0"}, {"dT", "0"}, {"min_score", "0"}});
    //             std::deque<int64_t> align_times = test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop", "argfile_last.override", false, true);
    //             std::ofstream fout_last("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop/align_time_last");
    //             fout_last << align_times.front() << '\n';
    //         }
    //         break;
    //     }
    //     std::ofstream fout_total("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_time_loop/total_times");
    //     fout_total << time(0) - time_tic << '\n';
    // }

    // override("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/zyj/zhengjin/4C_fold/run", "argfile", {{"align_threads_sz", "12"}, {"max_extract", "1"}, {"max_mega", "10000"}, {"diff_thres", "10"}, {"max_range", "10"}, {"max_round", "3"}, {"T", "-5"}, {"dT", "-5"}, {"min_score", "20"}});
    // test_Align("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/zyj/zhengjin/4C_fold/run", "argfile.override", false, true);

    // run_sim("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_single_cut", true);
    // run_sim("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/test_double_cut", true);

    std::deque<std::experimental::filesystem::path> dirs = {"/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Ct1-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Ct1-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Delta-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Delta-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Kappa-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Kappa-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Mu-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Mu-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Theta-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/Theta-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Ct1-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Ct1-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Delta-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Delta-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Kappa-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Kappa-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Mu-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Mu-Rep2",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Theta-Rep1",
    "/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Theta-Rep2"};

    int align_trheads_sz = 24, max_extract = 1, max_mega = 10000, diff_thres = 2, max_range = 10, min_seg_num = 1, max_seg_num = 3, block_size = 100, max_round = 4, min_score = 20;
    double T = -10, dT = -5;
    bool para = true;

    std::experimental::filesystem::path dir("/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/double/Ct1-Rep1");
    gradual_align(dir, "argfile", "argfile_last", dir.filename(), align_trheads_sz, max_extract, max_mega, diff_thres, max_range, min_seg_num, max_seg_num, block_size, max_round, T, dT, min_score, false, para);

    return 0;
}