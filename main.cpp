#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "headers/Track.h"
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <exception>
#include <unistd.h>
#include <ctime>
#include "headers/Help.h"

struct Parallel_Align
{
    int argc;
    char **argv;

    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    Graph graph;

    std::string read_file;
    std::deque<std::pair<std::string, std::string>> reads_total;
    int threads_sz = 24, max_extract = 1, diff_thres = 10, max_range = 10, min_seg_num = 0, max_seg_num = 0, block_size = 100;
    int64_t max_mega = 10000;
    int Omax;
    std::atomic<int> block_num;
    std::deque<std::string> mg_files;

    Parallel_Align(int argc_, char **argv_)
        : argc(argc_), argv(argv_), graph(argc, argv, file2seq, file2browheel)
    {
        for (int i = 1; i < argc; ++i)
        {
            if (!strcmp(argv[i], "---read"))
                read_file = argv[i+1];
            if (!strcmp(argv[i], "---threads_sz"))
                threads_sz = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---max_extract"))
                max_extract = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---diff_thres"))
                diff_thres = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---max_range"))
                max_range = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---min_seg_num"))
                min_seg_num = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---max_seg_num"))
                max_seg_num = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---max_mega"))
                max_mega = std::stoi(argv[i + 1]);
            if (!strcmp(argv[i], "---block_size"))
                block_size = std::stoi(argv[i + 1]);
        }

        Omax = 0;
        std::ifstream fin(read_file);
        int64_t preg = 0;
        while (fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
        {
            int64_t nowg = fin.tellg();
            int len = nowg - preg;
            if (len > Omax)
                Omax = len;
            preg = nowg;
        }
        --Omax;

        graph.draw("graph.gv");
    }

    void load_index()
    {
        for (auto &pair : file2browheel)
            pair.second.loadBroWheel(pair.first);
        for (EdgeGlobalCross &global : graph.global_crosses)
            if (file2browheel[global.name].reverse_complement != global.reverse_complement)
            {
                std::cerr << "the long reference " << global.name << " is required to be "
                          << global.reverse_complement << ", but the index is " << file2browheel[global.name].reverse_complement << '\n';
                exit(EXIT_FAILURE);
            }
        for (EdgeGlobalCircuit &global : graph.global_circuits)
            if (file2browheel[global.name].reverse_complement != global.reverse_complement)
            {
                std::cerr << "the long reference " << global.name << " is required to be "
                          << global.reverse_complement << ", but the index is " << file2browheel[global.name].reverse_complement << '\n';
                exit(EXIT_FAILURE);
            }
    }

    std::queue<std::pair<std::string, std::string>> fill_block(std::ifstream &fin)
    {
        std::queue<std::pair<std::string, std::string>> reads;
        while (reads.size() < block_size)
        {
            std::string name, read;
            if (std::getline(std::getline(fin, name), read))
            {
                reads.emplace(name, read);
                reads_total.emplace_back(name, read);
            }
            else
                break;
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
        Align align(argc, argv, file2seq, file2browheel, "mg", Omax);

        std::ifstream fin(read_file);
        for (std::queue<std::pair<std::string, std::string>> reads=fill_block(fin); !reads.empty(); reads=fill_block(fin))
            process_and_write_single(std::move(reads), align);
        fin.close();

        if (!align.Oname.empty())
            align.fout.write((char *)&align.max_id, sizeof(align.max_id));
        align.fout.close();
    }

    void parallel_align()
    {
        std::map<std::thread::id, Align> aligns;
        thread_pool threads1(threads_sz);
        auto thread_ids = threads1.get_ids();
        for (size_t i = 0; i < threads1.size(); ++i)
        {
            mg_files.emplace_back("mg" + std::to_string(i));
            aligns.emplace(std::piecewise_construct, std::forward_as_tuple(thread_ids[i]), std::forward_as_tuple(argc, argv, file2seq, file2browheel, mg_files.back(), Omax));
        }

        int max_block = 2 * threads_sz;
        block_num.store(0);
        std::deque<std::future<void>> futures;
        std::ifstream fin(read_file);
        for (std::queue<std::pair<std::string, std::string>> reads=fill_block(fin); !reads.empty(); reads=fill_block(fin))
        {
            while (block_num.load() >= max_block)
                ;
            futures.push_back(threads1.submit(std::bind(&Parallel_Align::process_and_write, this, std::move(reads), std::ref(aligns))));
            ++block_num;
        }
        fin.close();

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

    void track()
    {
        Track track(argc, argv, file2seq, file2browheel, max_extract, max_mega, diff_thres, max_range, min_seg_num, max_seg_num);
        track.ReadTrack(mg_files, reads_total);
    }
};

int main(int argc, char **argv)
{
    if (!strcmp(argv[1], "-h"))
    {
        help();
        return 0;
    }

    if (!strcmp(argv[1], "--index"))
    {
        BroWheel browheel;
        bool reverse_complement = false;
        std::string fasta;
        for (int i = 2; i < argc; ++i)
        {
            if (!strcmp(argv[i], "-r"))
                reverse_complement = true;
            if (!strcmp(argv[i], "-f"))
                fasta = argv[i + 1];
        }
        browheel.readin(fasta, true, reverse_complement);
        thread_pool threads1(3), thread2(1);
        browheel.index(150000000, threads1, thread2);
        browheel.saveBroWheel();
        return 0;
    }

    Parallel_Align parallel_align(argc, argv);
    parallel_align.load_index();
    parallel_align.parallel_align();
    parallel_align.track();

    return 0;
}
