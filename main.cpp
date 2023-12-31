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

QUERYSIZE set_Omax(std::string input)
{
    QUERYSIZE Omax = 0;
    for (std::ifstream fin(input); fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n').good();)
    {
        SIZETYPE preg = fin.tellg();
        QUERYSIZE len = fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n').tellg() - preg - 1;
        if (len > Omax)
            Omax = len;
    }

    return Omax;
}

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
        ("threads_sz", boost::program_options::value<SIZETYPE>()->default_value(24), "# of align threads");

    track_options.add_options()
        ("max_extract", boost::program_options::value<SIZETYPE>()->default_value(1), "max best maps to extract")
        ("max_mega", boost::program_options::value<SIZETYPE>()->default_value(10000), "max steps to track map graph")
        ("diff_thres", boost::program_options::value<SIZETYPE>()->default_value(10), "min bps to distinct two maps")
        ("max_range", boost::program_options::value<SIZETYPE>()->default_value(10), "max map pos to show per seg")
        ("min_seg_num", boost::program_options::value<SIZETYPE>()->default_value(0), "min map segments, 0 means no restriction")
        ("max_seg_num", boost::program_options::value<SIZETYPE>()->default_value(0), "max map segments, 0 means no restriction");

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
        SIZETYPE revref_sz = std::filesystem::file_size(RevRefFile);
        NUCTYPE *revref = new NUCTYPE[revref_sz];
        std::ifstream fin(RevRefFile);
        fin.read((char*)revref, revref_sz);
        gsufsortSA(RevRefFile+".sa", revref, revref_sz);
        RankVec bwtRank;
        bwtRank.bwt_sz = revref_sz;
        bwtRank.bwt = SA2bwt(RevRefFile+".sa", 10000, revref);
        delete[] revref;
        bwtRank.count();
        bwtRank.saveRankVec(RevRefFile+".bwt", RevRefFile+".rnk");

        // std::cerr << "check bwt:" << check_bwt(bwtRank, revref) << '\n'; // debug

        return EXIT_SUCCESS;
    }

    std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> file2short;
    if (vm.count("shorts"))
        for (std::string file : vm["shorts"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2short.count(file))
                continue;
            file2short[file].second = std::filesystem::file_size(file);
            file2short[file].first.reset(new NUCTYPE[file2short[file].second]);
            std::ifstream fin(file);
            fin.read((char*)file2short[file].first.get(), file2short[file].second);
        }
    
    std::map<std::string, RankVec> file2rankvec;
    std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SIZETYPE>> file2long;
    if (vm.count("longs"))
        for (std::string file : vm["longs"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2rankvec.count(file))
                continue;
            file2rankvec[file].loadRankVec(file+".bwt", file+".rnk");
        }

    {
        graph_t graph;
        construct_graph(graph, vm, file2short);
        draw(graph, file2short);
    }

    QUERYSIZE Omax = set_Omax(vm["input"].as<std::string>());

    std::unique_ptr<char []> inbuf(new char[10*1024*1024]);
    std::ifstream fin(vm["input"].as<std::string>());
    fin.rdbuf()->pubsetbuf(inbuf.get(), 10*1024*1024);
    std::mutex mtx;
    std::deque<std::thread> threads;
    std::deque<Align> aligns;
    std::deque<std::string> mg_files;
    for (SIZETYPE i = 0; i < vm["threads_sz"].as<SIZETYPE>(); ++i)
    {
        mg_files.emplace_back("mg" + std::to_string(i));
        aligns.emplace_back(mtx, fin, vm, file2short, file2rankvec, mg_files.back(), Omax);
        threads.emplace_back(&Align::run, &aligns.back());
    }
    for (SIZETYPE i = 0; i < vm["threads_sz"].as<SIZETYPE>(); ++i)
        threads[i].join();

    // Align align(mtx, fin, vm, file2short, file2rankvec, "mg", Omax); // debug
    // align.run(); // debug

    // std::deque<std::string> mg_files; // debug
    // for (SIZETYPE i = 0; i < vm["threads_sz"].as<SIZETYPE>(); ++i) // debug
    //     mg_files.emplace_back("mg" + std::to_string(i)); // debug

    Track track(vm, file2short, file2rankvec); 
    track.ReadTrack(mg_files);


    return EXIT_SUCCESS;
}
