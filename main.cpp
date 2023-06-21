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
#include <boost/program_options.hpp>

struct Parallel_Align
{
    int argc;
    char **argv;

    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    Graph graph;

    std::string input;
    int threads_sz, max_extract, diff_thres, max_range, min_seg_num, max_seg_num, block_size;
    int64_t max_mega;
    int Omax;
    std::atomic<int> block_num;
    std::deque<std::string> mg_files;

    Parallel_Align(int argc_, char **argv_)
        : argc(argc_), argv(argv_), graph(argc, argv, file2seq, file2browheel)
    {
    }

    void set_Omax_graph_draw()
    {
        Omax = 0;
        std::ifstream fin(input);
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
    }

    std::queue<std::pair<std::string, std::vector<int8_t>>> fill_block(std::ifstream &fin)
    {
        std::queue<std::pair<std::string, std::vector<int8_t>>> reads;
        while (reads.size() < block_size)
        {
            std::string name, read;
            if (std::getline(std::getline(fin, name), read))
            {
                reads.emplace(name, std::vector<int8_t>());
                for (int i = 0; i < read.size(); ++i)
                    reads.back().second.push_back(BroWheel::base2int[read[i]]);
            }
            else
                break;
        }
        return reads;
    }

    void process_and_write_single(std::queue<std::pair<std::string, std::vector<int8_t>>> &&reads, Align &align)
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

        std::ifstream fin(input);
        for (std::queue<std::pair<std::string, std::vector<int8_t>>> reads=fill_block(fin); !reads.empty(); reads=fill_block(fin))
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
        std::ifstream fin(input);
        for (std::queue<std::pair<std::string, std::vector<int8_t>>> reads=fill_block(fin); !reads.empty(); reads=fill_block(fin))
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

    void process_and_write(std::queue<std::pair<std::string, std::vector<int8_t>>> reads, std::map<std::thread::id, Align> &aligns)
    {
        process_and_write_single(std::move(reads), aligns.at(std::this_thread::get_id()));
        --block_num;
    }

    void track()
    {
        Track track(argc, argv, file2seq, file2browheel, max_extract, max_mega, diff_thres, max_range, min_seg_num, max_seg_num);
        track.ReadTrack(mg_files, input);
    }
};

int main(int argc, char **argv)
{
    Parallel_Align parallel_align(argc, argv);

    boost::program_options::options_description general_options("general options"), index_options("index options"), parallel_options("parallel options"), output_options("output options"), map_options("map options"), all_options("all options");
    general_options.add_options()
        ("help,h", "help screen")
        ("input", boost::program_options::value<std::string>(&parallel_align.input), "input");
    index_options.add_options()
        ("index", boost::program_options::value<std::string>(), "reference file to index");
    parallel_options.add_options()
        ("threads_sz", boost::program_options::value<int>(&parallel_align.threads_sz)->default_value(24), "# of align threads")
        ("block_size", boost::program_options::value<int>(&parallel_align.block_size)->default_value(10), "size of reads batch");
    output_options.add_options()
        ("max_range", boost::program_options::value<int>(&parallel_align.max_range)->default_value(10), "max map pos to show per seg")
        ("max_extract", boost::program_options::value<int>(&parallel_align.max_extract)->default_value(1), "max best maps to extract")
        ("max_mega", boost::program_options::value<int64_t>(&parallel_align.max_mega)->default_value(10000), "max steps to track map graph");
    map_options.add_options()
        ("diff_thres", boost::program_options::value<int>(&parallel_align.diff_thres)->default_value(10), "min bps to distinct two maps")
        ("min_seg_num", boost::program_options::value<int>(&parallel_align.min_seg_num)->default_value(0), "min map segments, 0 means no restriction")
        ("max_seg_num", boost::program_options::value<int>(&parallel_align.max_seg_num)->default_value(0), "max map segments, 0 means no restriction");

    all_options.add(general_options).add(index_options).add(parallel_options).add(output_options).add(map_options);
    boost::program_options::variables_map vm;
    boost::program_options::command_line_parser parser{argc, argv};
    boost::program_options::store(parser.options(all_options).allow_unregistered().run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("help"))
    {
        std::cout << all_options << '\n';
        return EXIT_SUCCESS;
    }

    if (vm.count("index"))
    {
        BroWheel browheel;
        browheel.readin(vm["index"].as<std::string>(), true);
        thread_pool threads1(3), thread2(1);
        browheel.index(150000000, threads1, thread2);
        // browheel.index(50, threads1, thread2);
        browheel.saveBroWheel();
        return EXIT_SUCCESS;
    }

    parallel_align.set_Omax_graph_draw();
    parallel_align.load_index();
    parallel_align.parallel_align();
    parallel_align.track();

    // BroWheel browheel;
    // browheel.loadBroWheel("test.fa");
    // std::cerr << browheel.bwt.size() << '\n';
    // std::ofstream fout("test.fa.my.bwt");
    // for (int i=0; i<browheel.bwt.size(); ++i)
    //     fout << BroWheel::int2base[browheel.bwt[i]];
    // fout << '\n';
    // fout.close();

    return EXIT_SUCCESS;
}
