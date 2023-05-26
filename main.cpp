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

void help()
{
    std::cout << "--read_files file1 file2 ...\n"
    << "--nodes node_name1 gap_open_penalty(GOP)1 gap_extend_penalty(GEP)1 node_name2...\n"
    << "--roots root_name1 root_name2 ...\n"
    << "--targets target_name1 target_name2 ...\n"
    << "--locals <default_gamma | 49_pair_scores1> local_file1 GOP_in_ref1 GEP_in_ref1 GOP_in_query1 GEP_in_query1 basic_penalty1 GOP_in_ref_begin1 GEP_in_ref_begin1 GOP_in_ref_end1 GEP_in_ref_end1 tail1 head1 <default_gamma | 49_pair_scores2> local_file2...\n"
    << "--globals <default_gamma | 49_pair_scores1> global_file1 GOP_in_ref1 GEP_in_ref1 GOP_in_query1 GEP_in_query1 basic_penalty1 tail1 head1 <default_gamma | 49_pair_scores2> global_file2...\n";
}

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
    auto &align = aligns.at(std::this_thread::get_id());
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

std::deque<std::future<void>> load_blocks(std::deque<std::string> read_files, thread_pool &threads1, std::map<std::thread::id, Align> &aligns, size_t max_block, size_t block_size)
{
    std::string name;
    std::string read;
    auto iter = read_files.begin();
    std::ifstream fin(*iter);
    std::deque<std::future<void>> futures;
    while (iter != read_files.end())
    {
        std::queue<std::pair<std::string, std::string>> reads;
        while (reads.size() < block_size && iter != read_files.end())
        {
            if (std::getline(std::getline(fin, name), read))
            {
                reads.emplace(name, read);
                if (name[0]=='@')
                    fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n').ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
            else
            {
                fin.close();
                fin.open(*(++iter));
            }
        }
        std::unique_lock lck(mtx);
        cv.wait(lck, [max_block]
                { return block_num < max_block; });
        futures.push_back(threads1.submit(std::bind(process_and_write, std::move(reads), std::ref(aligns))));
        ++block_num;
    }
    fin.close();
    return futures;
}

void load_blocks(std::deque<std::string> read_files, Align &align, size_t block_size)
{
    std::string name;
    std::string read;
    auto iter = read_files.begin();
    std::ifstream fin(*iter);
    while (iter != read_files.end())
    {
        std::queue<std::pair<std::string, std::string>> reads;
        while (reads.size() < block_size && iter != read_files.end())
        {
            if (std::getline(std::getline(fin, name), read))
            {
                reads.emplace(name, read);
                if (name[0]=='@')
                    fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n').ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
            else
            {
                fin.close();
                fin.open(*(++iter));
            }
        }
        ++block_num;
        process_and_write_single(std::move(reads), align);
    }
    fin.close();
}

std::tuple<int, char **> read_argfile(std::string argfile, std::vector<std::string> &args)
{
    std::string tmp;
    std::ifstream fin(argfile);
    while (fin >> tmp)
    {
        if (tmp.back()==',')
            args.push_back(tmp.substr(1, tmp.size() - 3));
        else
            args.push_back(tmp.substr(1, tmp.size() - 2));
    }
    fin.close();
    int argc = args.size() + 1;
    char **argv = new char *[argc];
    for (int i = 1; i < argc; ++i)
        argv[i] = (char *)args[i - 1].c_str();

    return std::make_tuple(argc, argv);
}

std::deque<std::string> load_index(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel)
{
    Graph graph(argc, argv, file2seq, file2browheel);
    graph.draw("graph.gv");
    for (auto &pair : file2seq)
        pair.second.readin(pair.first);
    for (auto &pair : file2browheel)
        pair.second.loadBroWheel(pair.first);
    std::deque<std::string> read_files;
    for (int i = 1; i < argc; ++i)
        if (!strcmp(argv[i], "--read_files"))
        {
            while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                read_files.push_back(argv[i]);
            --i;
        }

    return read_files;
}

std::deque<std::string> read_reference_and_index(int argc, char **argv, int ll, thread_pool &threads1, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel, bool reverse_complement)
{
    Graph graph(argc, argv, file2seq, file2browheel);
    graph.draw("graph.gv");
    for (auto &pair : file2seq)
        pair.second.readin(pair.first);
    for (auto &pair : file2browheel)
        pair.second.readin(pair.first, true, reverse_complement);
    thread_pool thread2(1);
    for (auto &pair : file2browheel)
    {
        pair.second.index(ll, threads1, thread2);
        pair.second.saveBroWheel();
    }
    std::deque<std::string> read_files;
    for (int i = 1; i < argc; ++i)
        if (!strcmp(argv[i], "--read_files"))
        {
            while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                read_files.push_back(argv[i]);
            --i;
        }

    return read_files;
}

std::deque<std::string> parallel_align(size_t chunk_sz, thread_pool &threads1, int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel, std::deque<std::string> read_files, std::string run_name)
{
    std::map<std::thread::id, Align> aligns;
    auto thread_ids = threads1.get_ids();
    std::deque<std::string> mg_files;
    for (size_t i = 0; i < threads1.size(); ++i)
    {
        mg_files.emplace_back(run_name + std::to_string(i) + ".mg");
        aligns.emplace(std::piecewise_construct, std::forward_as_tuple(thread_ids[i]), std::forward_as_tuple(argc, argv, file2seq, file2browheel, chunk_sz, mg_files.back()));
    }

    size_t block_size = 100, max_block = 2 * threads1.size();
    std::deque<std::future<void>> futures = load_blocks(read_files, threads1, aligns, max_block, block_size); // mt
    for (auto &future : futures)                                                                             // mt
        future.wait();                                                                                       // mt
    // load_blocks(read_files, aligns.begin()->second, block_size); // st

    for (auto &pair : aligns)
    {
        auto &align = pair.second;
        if (!align.Oname.empty())
            align.fout.write((char *)&align.max_id, sizeof(align.max_id));
        align.fout.close();
    }

    return mg_files;
}

std::deque<std::string> get_dirs()
{   
    const std::experimental::filesystem::path current_path(".");
    std::experimental::filesystem::directory_iterator end_itr; // default construction yields past-the-end
    std::deque<std::string> dirs;
    for ( std::experimental::filesystem::directory_iterator itr( current_path ); itr != end_itr; ++itr )
        if ( std::experimental::filesystem::is_directory(itr->status()) )
            dirs.push_back(itr->path().filename().string());
    return dirs;
}

std::deque<std::string> get_argfiles()
{   
    const std::experimental::filesystem::path current_path(".");
    std::experimental::filesystem::directory_iterator end_itr; // default construction yields past-the-end
    std::deque<std::string> argfiles;
    for ( std::experimental::filesystem::directory_iterator itr( current_path ); itr != end_itr; ++itr )
        if ( ! std::experimental::filesystem::is_directory(itr->status()) )
        {
            std::string str_tmp=itr->path().filename().string();
            if (str_tmp.size()>8 && str_tmp.substr(str_tmp.size()-8,8)==".argfile")
                argfiles.push_back(str_tmp);
        }
    return argfiles;
}

void run_sim(std::string sim_dir);

void test_RandomReads(std::string argfile, int threads1_sz);
void index_genome(std::string argfile, int threads1_sz);
void test_CRISPR(std::string argfile, int threads1_sz);
void test_single_double(std::string argfile, int threads1_sz);

int main(int argc, char **argv)
{
    // test_RandomReads("argfile", 24);

    chdir("./test_CRISPR_cross");
    // index_genome("argfile", 3);
    test_CRISPR("argfile", 3);

    // run_sim("./test_single_cut");
    // run_sim("./test_double_cut");

    return 0;
}

void run_sim(std::string sim_dir)
{
    chdir(sim_dir.c_str());
    std::deque<std::string> dirs=get_dirs();
    for (auto dir : dirs)
    {
        std::string tmp="./"+dir;
        chdir(tmp.c_str());
        std::deque<std::string> argfiles=get_argfiles();
        for (auto &argfile : argfiles)
            test_single_double(argfile, 12);
        chdir("..");
    }
    chdir("..");
}

void test_RandomReads(std::string argfile, int threads1_sz)
{
    int n_sz = 5, r_sz = 2, t_sz = 1, max_e_sz = 8, lseqlb = 100, lsequb = 200, gseqlb = 10000, gsequb = 20000, aseqlb = 50, asequb = 100, seq_num = 10000, head_in = 10, tail_in = 10;
    bool acyclic = false;
    double gpro = 0.5, rpro = 0, nve = 0, nue = 0, ve = -5, ue = -2, vf = -5, uf = -2, T = -10, vfp = 0, ufp = 0, vfm = 0, ufm = 0, mat = 1, mis = -3, indel_rate = 0.005, mut_rate = 0.005, apro = 0.5;

    random_DG(n_sz, r_sz, t_sz, max_e_sz, "argfile", acyclic, gpro, rpro, nve, nue, ve, ue, vf, uf, T, vfp, ufp, vfm, ufm, mat, mis, "local_file", "global_file", lseqlb, lsequb, gseqlb, gsequb, aseqlb, asequb, seq_num, "read_file", "truth_file", indel_rate, mut_rate, head_in, tail_in, apro);

    int argc;
    char **argv;
    std::vector<std::string> args;
    std::tie(argc, argv) = read_argfile(argfile, args);

    thread_pool threads1(threads1_sz);
    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    std::deque<std::string> read_files = read_reference_and_index(argc, argv, 300, threads1, file2seq, file2browheel, false);

    std::deque<std::string> mg_files = parallel_align(128 * 1024 * 1024, threads1, argc, argv, file2seq, file2browheel, read_files, "Random");

    Track track(argc, argv, file2seq, file2browheel, 128 * 1024 * 1024);
    track.ReadTrack(mg_files, read_files, "Random", INT64_MAX, 1);
    track.extract("Random", 1);
}

void index_genome(std::string argfile, int threads1_sz)
{
    int argc;
    char **argv;
    std::vector<std::string> args;
    std::tie(argc, argv) = read_argfile(argfile, args);

    thread_pool threads1(threads1_sz);
    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    std::deque<std::string> read_files = read_reference_and_index(argc, argv, 150000000, threads1, file2seq, file2browheel, true);
    // std::deque<std::string> read_files = read_reference_and_index(argc, argv, 150000000, threads1_sz, file2seq, file2browheel, true);
}

void test_CRISPR(std::string argfile, int threads1_sz)
{
    int argc;
    char **argv;
    std::vector<std::string> args;
    std::tie(argc, argv) = read_argfile(argfile, args);

    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    std::deque<std::string> read_files = load_index(argc, argv, file2seq, file2browheel);

    thread_pool threads1(threads1_sz);
    std::deque<std::string> mg_files = parallel_align(128 * 1024 * 1024, threads1, argc, argv, file2seq, file2browheel, read_files, "CRISPR");

    Track track(argc, argv, file2seq, file2browheel, 128 * 1024 * 1024);
    track.ReadTrack(mg_files, read_files, "CRISPR", INT64_MAX, 1);
    track.extract("CRISPR", 1);
}

void test_single_double(std::string argfile, int threads1_sz)
{
    int argc;
    char **argv;
    std::vector<std::string> args;
    std::tie(argc, argv) = read_argfile(argfile, args);

    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    std::deque<std::string> read_files = load_index(argc, argv, file2seq, file2browheel);

    thread_pool threads1(threads1_sz);
    std::deque<std::string> mg_files = parallel_align(128 * 1024 * 1024, threads1, argc, argv, file2seq, file2browheel, read_files, "sim");

    Track track(argc, argv, file2seq, file2browheel, 128 * 1024 * 1024);
    track.ReadTrack(mg_files, read_files, "sim", INT64_MAX, 10);
    track.extract("sim", 10);
}