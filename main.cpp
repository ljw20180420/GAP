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


int64_t set_Omax(std::string input)
{
    int64_t Omax = 0;
    std::ifstream fin(input);
    int64_t preg = 0;
    while (fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
    {
        int64_t nowg = fin.tellg(); // tellg returns -1 for exception, so nowg cannot be uint64_t
        int64_t len = nowg - preg;
        if (len > Omax)
            Omax = len;
        preg = nowg;
    }
    --Omax; // empty input makes Omax=-1, so Omax cannot be uint64_t

    return Omax;
}


struct Parallel_Align
{
    boost::program_options::variables_map &vm;
    std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2short;
    std::map<std::string, RankVec> &file2rankvec;
    std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2long;
    std::atomic<uint64_t> block_num;

    Parallel_Align(boost::program_options::variables_map &vm_, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2short_, std::map<std::string, RankVec> &file2rankvec_, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2long_) : vm(vm_), file2short(file2short_), file2rankvec(file2rankvec_), file2long(file2long_)
    {}

    std::queue<std::pair<std::string, std::vector<uint8_t>>> fill_block(std::ifstream &fin)
    {
        std::queue<std::pair<std::string, std::vector<uint8_t>>> reads;
        while (reads.size() < vm["block_size"].as<uint64_t>())
        {
            std::string name, read;
            if (std::getline(std::getline(fin, name), read))
            {
                reads.emplace(name, std::vector<uint8_t>());
                for (uint64_t i = 0; i < read.size(); ++i)
                    reads.back().second.push_back(base2int[read[i]]);
            }
            else
                break;
        }
        return reads;
    }

    void process_and_write_single(std::queue<std::pair<std::string, std::vector<uint8_t>>> &&reads, Align &align)
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

    std::deque<std::string> single_align(int64_t Omax)
    {
        Align align(vm, file2short, file2rankvec, file2long, "mg", Omax);

        std::ifstream fin(vm["input"].as<std::string>());
        for (std::queue<std::pair<std::string, std::vector<uint8_t>>> reads=fill_block(fin); !reads.empty(); reads=fill_block(fin))
            process_and_write_single(std::move(reads), align);
        fin.close();

        if (!align.Oname.empty())
            align.fout.write((char *)&align.max_id, sizeof(align.max_id));
        align.fout.close();

        return std::deque<std::string>{"mg"};
    }

    std::deque<std::string> parallel_align(int64_t Omax)
    {
        std::map<std::thread::id, Align> aligns;
        thread_pool threads1(vm["threads_sz"].as<uint64_t>());
        auto thread_ids = threads1.get_ids();
        std::deque<std::string> mg_files;
        for (uint64_t i = 0; i < threads1.size(); ++i)
        {
            mg_files.emplace_back("mg" + std::to_string(i));
            aligns.emplace(std::piecewise_construct, std::forward_as_tuple(thread_ids[i]), std::forward_as_tuple(vm, file2short, file2rankvec, file2long, mg_files.back(), Omax));
        }

        uint64_t max_block = 2 * vm["threads_sz"].as<uint64_t>();
        block_num.store(0);
        std::deque<std::future<void>> futures;
        std::ifstream fin(vm["input"].as<std::string>());
        for (std::queue<std::pair<std::string, std::vector<uint8_t>>> reads=fill_block(fin); !reads.empty(); reads=fill_block(fin))
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
        return mg_files;
    }

    void process_and_write(std::queue<std::pair<std::string, std::vector<uint8_t>>> reads, std::map<std::thread::id, Align> &aligns)
    {
        process_and_write_single(std::move(reads), aligns.at(std::this_thread::get_id()));
        --block_num;
    }
};

int main(int argc, char **argv)
{
    boost::program_options::options_description general_options("general options"), graph_options("graph options"), align_options("align options"), track_options("track options");

    general_options.add_options()
        ("help,h", "help screen")
        ("index", boost::program_options::value<std::string>(), "reference file to index");

    graph_options.add_options()
        ("nodes", boost::program_options::value<std::vector<std::string>>()->multitoken()->
        zero_tokens()->composing(), "graph nodes: name,is_root,is_target,v,u")
        ("shorts", boost::program_options::value<std::vector<std::string>>()->multitoken()->
        zero_tokens()->composing(), "graph short edges: name,tail,head,match,mismatch,v,u,T,min_score,vb,ub")
        ("longs", boost::program_options::value<std::vector<std::string>>()->multitoken()->
        zero_tokens()->composing(), "graph long edges: name,tail,head,match,mismatch,v,u,T,min_score");

    align_options.add_options()
        ("input", boost::program_options::value<std::string>(), "input")
        ("threads_sz", boost::program_options::value<uint64_t>()->default_value(24), "# of align threads")
        ("block_size", boost::program_options::value<uint64_t>()->default_value(10), "size of reads batch");

    track_options.add_options()
        ("max_extract", boost::program_options::value<int>()->default_value(1), "max best maps to extract")
        ("max_mega", boost::program_options::value<int64_t>()->default_value(10000), "max steps to track map graph")
        ("diff_thres", boost::program_options::value<int>()->default_value(10), "min bps to distinct two maps")
        ("max_range", boost::program_options::value<int>()->default_value(10), "max map pos to show per seg")
        ("min_seg_num", boost::program_options::value<int>()->default_value(0), "min map segments, 0 means no restriction")
        ("max_seg_num", boost::program_options::value<int>()->default_value(0), "max map segments, 0 means no restriction");

    general_options.add(graph_options).add(align_options).add(track_options);
    boost::program_options::variables_map vm;
    boost::program_options::command_line_parser parser{argc, argv};
    boost::program_options::store(parser.options(general_options).allow_unregistered().run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("help"))
    {
        std::cout << general_options << '\n';
        return EXIT_SUCCESS;
    }

    if (vm.count("index"))
    {
        std::string RevRefFile = vm["index"].as<std::string>();
        uint64_t revref_sz = std::filesystem::file_size(RevRefFile);
        uint8_t *revref = (uint8_t *) malloc(revref_sz);
        std::ifstream fin(RevRefFile);
        fin.read((char*)revref, revref_sz);
        gsufsortSA(RevRefFile+".sa", revref, revref_sz);
        RankVec bwtRank;
        bwtRank.bwt_sz = revref_sz;
        bwtRank.bwt = SA2bwt(RevRefFile+".sa", 10000, revref);
        std::free(revref);
        bwtRank.count();
        bwtRank.saveRankVec(RevRefFile+".bwt", RevRefFile+".rnk");
        return EXIT_SUCCESS;
    }

    std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> file2short;
    if (vm.count("shorts"))
        for (std::string file : vm["shorts"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2short.count(file))
                continue;
            file2short[file].second = std::filesystem::file_size(file);
            file2short[file].first.reset((uint8_t *) malloc(file2short[file].second));
            std::ifstream fin(file);
            fin.read((char*)file2short[file].first.get(), file2short[file].second);
        }
    
    std::map<std::string, RankVec> file2rankvec;
    std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> file2long;
    if (vm.count("longs"))
        for (std::string file : vm["longs"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2rankvec.count(file))
                continue;
            file2rankvec[file].loadRankVec(file+".bwt", file+".rnk");

            file2long[file].second = std::filesystem::file_size(file);
            file2long[file].first.reset((uint8_t *) malloc(file2long[file].second));
            std::ifstream fin(file);
            fin.read((char *)file2long[file].first.get(), file2long[file].second);
        }

    {
        Graph graph(vm, file2short, file2rankvec, file2long);
        graph.draw("graph.gv");
    }

    int64_t Omax = set_Omax(vm["input"].as<std::string>());

    Parallel_Align parallel_align(vm, file2short, file2rankvec, file2long);
    std::deque<std::string> mg_files = parallel_align.parallel_align(Omax);
    // std::deque<std::string> mg_files = parallel_align.single_align(Omax);

    // std::deque<std::string> mg_files; // debug
    // for (uint64_t i = 0; i < vm["threads_sz"].as<uint64_t>(); ++i) // debug
    //     mg_files.emplace_back("mg" + std::to_string(i)); // debug

    Track track(vm, file2short, file2rankvec, file2long); 
    track.ReadTrack(mg_files);

    return EXIT_SUCCESS;
}
