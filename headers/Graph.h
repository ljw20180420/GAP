#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <cfloat>
#include <list>
#include "BroWheel.h"
#include <boost/program_options.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>

using DOTTYPE = int16_t;
using QUERYSIZE = uint32_t;
using SHORTSIZE = uint32_t;
using IDTYPE = uint64_t;
using SCORETYPE = int32_t;
using SOURCESIZE = uint32_t;
using BITSTYPE = uint8_t; // this is used for track the alignment, e.g. 8bits (uint8_t) can track 8 maxmial sources
constexpr static const SCORETYPE inf = std::numeric_limits<SCORETYPE>::max() / 2;

struct SNC
{
    SCORETYPE E;
    SCORETYPE F0;
    SCORETYPE G0;
    SCORETYPE hatG;
    NUCTYPE c;
    NUCTYPE idl;
    QUERYSIZE lambda;
    SIZETYPE sr1;
    SIZETYPE sr2;
    IDTYPE cid; // 0 means not try to add child yet
    SNC *itp;
    SNC *jump;

    SNC(SCORETYPE E_, SCORETYPE F0_, SCORETYPE G0_, SCORETYPE hatG_, NUCTYPE c_, NUCTYPE idl_, QUERYSIZE lambda_, SIZETYPE sr1_, SIZETYPE sr2_, SIZETYPE cid_, SNC *itp_, SNC *jump_)
        : E(E_), F0(F0_), G0(G0_), hatG(hatG_), c(c_), idl(idl_), lambda(lambda_), sr1(sr1_), sr2(sr2_), cid(cid_), itp(itp_), jump(jump_)
    {
    }
};

struct Dot
{
    static const DOTTYPE DotQ = -2, DotAbar = -1;

    DOTTYPE n; // n determines the type of Dot
    SIZETYPE s;
    QUERYSIZE w;
    SCORETYPE val;
    SIZETYPE fs; // first source
    SOURCESIZE s_sz;
    QUERYSIZE lambda;
    IDTYPE id;

    static DOTTYPE DOTTYPE2NODEIDX(DOTTYPE nidx)
    {
        return (-nidx + DotQ - 1) / 2;
    }

    static DOTTYPE NODEIDX2DOTTYPEB(DOTTYPE nidx)
    {
        return -2 * nidx + DotQ - 1;
    }

    static DOTTYPE NODEIDX2DOTTYPEA(DOTTYPE nidx)
    {
        return -2 * nidx + DotQ - 2;
    }
};

struct TrackTree
{
    DOTTYPE n;
    QUERYSIZE W, lambda;
    std::deque<Dot> dots;
    int64_t idx, pidx;

    struct TrackNode
    {
        int64_t cidxs[5] = {-1, -1, -1, -1, -1};
        QUERYSIZE tau = 0;
    };

    std::deque<TrackNode> tracknodes;

    void initialize(DOTTYPE n_, QUERYSIZE W_)
    {
        lambda = 0;
        n = n_;
        W = W_;
        tracknodes.clear();
        emplace_back(0);
    }

    void emplace_back(SIZETYPE s)
    {
        tracknodes.emplace_back();
        dots.resize(tracknodes.size() * 3 * (W + 1));
        SIZETYPE k = (tracknodes.size() - 1) * 3 * (W + 1);
        for (SIZETYPE i = 0; i < 3; ++i)
            for (SIZETYPE w = 0; w <= W; ++w, ++k)
            {
                dots[k].n = n;
                dots[k].s = s;
                dots[k].w = w;
                dots[k].lambda = lambda;
            }
    }

    void goto_child(NUCTYPE c, SIZETYPE s)
    {
        ++lambda;
        pidx = idx;
        if (tracknodes[pidx].cidxs[c - 2] < 0)
        {
            idx = tracknodes.size();
            tracknodes[pidx].cidxs[c - 2] = idx;
            emplace_back(s);
        }
        else
            idx = tracknodes[pidx].cidxs[c - 2];
    }

    Dot & dotE(SIZETYPE w)
    {
        return dots[3 * idx * (W + 1) + w];
    }

    Dot & dotF(SIZETYPE w)
    {
        return dots[(3 * idx + 1) * (W + 1) + w];
    }

    Dot & dotG(SIZETYPE w)
    {
        return dots[(3 * idx + 2) * (W + 1) + w];
    }

    Dot & dotPE(SIZETYPE w)
    {
        return dots[3 * pidx * (W + 1) + w];
    }

    Dot & dotPF(SIZETYPE w)
    {
        return dots[(3 * pidx + 1) * (W + 1) + w];
    }

    Dot & dotPG(SIZETYPE w)
    {
        return dots[(3 * pidx + 2) * (W + 1) + w];
    }
};

struct Node;

struct Edge
{
    // basic information
    DOTTYPE n;
    std::string name;

    // basic score
    SCORETYPE gamma[RankVec::sigma][RankVec::sigma] = {
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -3, -3, -3, -3, -3},
        {-inf, -inf, -3, 1, -3, -3, -3},
        {-inf, -inf, -3, -3, 1, -3, -3},
        {-inf, -inf, -3, -3, -3, 1, -3},
        {-inf, -inf, -3, -3, -3, -3, 1}};
    SCORETYPE ve, ue, vf, uf, T, min_score;
    
    // short specific score
    SCORETYPE vfp, ufp, vfm, ufm, gfm;
    std::vector<SCORETYPE> gf, gfpT;

    // cross short specific data
    SCORETYPE **Fvals = NULL;

    // short specific data
    SCORETYPE **Evals = NULL, **Gvals = NULL;

    // circuit short specific data
    SCORETYPE **F0vals = NULL, **G0vals = NULL, **DXvals = NULL;
    SCORETYPE *Dvals;

    // circuit specific data
    SCORETYPE **D0vals = NULL;

    // circuit long specific data
    std::deque<SNC> sncs;

    // long specific data for track
    TrackTree tracktree;
    SCORETYPE *tailAvals;

    ~Edge()
    {
        if (Fvals)
            delete[] Fvals;
        if (Evals)
            delete[] Evals;
        if (Gvals)
            delete[] Gvals;
        if (F0vals)
            delete[] F0vals;
        if (G0vals)
            delete[] G0vals;
        if (DXvals)
            delete[] DXvals;
        if (D0vals)
            delete[] D0vals;
        if (tailAvals)
            delete[] tailAvals;
    }
};

namespace boost {
    enum vertex_Node_t {vertex_Node = 111};
    BOOST_INSTALL_PROPERTY(vertex, Node);
    enum vertex_comp_t {vertex_comp = 112};
    BOOST_INSTALL_PROPERTY(vertex, comp);
    enum vertex_compsz_t {vertex_compsz = 113};
    BOOST_INSTALL_PROPERTY(vertex, compsz);
    enum edge_Edge_t {edge_Edge = 114};
    BOOST_INSTALL_PROPERTY(edge, Edge);
}
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_compsz_t, SIZETYPE, boost::property<boost::vertex_comp_t, SIZETYPE, boost::property<boost::vertex_Node_t, Node>>>, boost::property<boost::edge_Edge_t, Edge>> graph_t;

struct Node
{
    bool is_root, is_target;
    std::string name;
    SCORETYPE ve, ue;

    SCORETYPE ** Avals = NULL;
    SCORETYPE *Bvals, *pAbarval;

    std::deque<SCORETYPE *>* AdeltaDot = NULL;

    struct GlobalSuffix
    {
        boost::graph_traits<graph_t>::edge_iterator ei;
        SIZETYPE start;
        QUERYSIZE lambda;

        GlobalSuffix(boost::graph_traits<graph_t>::edge_iterator ei_, SIZETYPE start_, QUERYSIZE lambda_)
            : ei(ei_), start(start_), lambda(lambda_)
        {
        }
    };

    std::deque<GlobalSuffix>* AdeltaGlobal = NULL;

    ~Node()
    {
        if (Avals)
            delete[] Avals;
        if (AdeltaDot)
            delete[] AdeltaDot;
        if (AdeltaGlobal)
            delete[] AdeltaGlobal;
    }
};

void get_affine(std::vector<SCORETYPE> &g, SHORTSIZE seqlen, SCORETYPE initial, SCORETYPE u, SCORETYPE v)
{
    g.resize(seqlen + 1);
    g[0] = initial;
    for (SIZETYPE s = 0; s < seqlen; ++s)
        if (s == 0)
            g[s + 1] = g[s] + v;
        else
            g[s + 1] = g[s] + u;
}

SIZETYPE construct_graph(graph_t &graph, boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short)
{
    boost::property_map<graph_t, boost::vertex_Node_t>::type node_map = boost::get(boost::vertex_Node, graph);
    boost::property_map<graph_t, boost::edge_Edge_t>::type edge_map = boost::get(boost::edge_Edge, graph);
    boost::property_map<graph_t, boost::vertex_comp_t>::type comp_map = boost::get(boost::vertex_comp, graph);
    boost::property_map<graph_t, boost::vertex_compsz_t>::type compsz_map = boost::get(boost::vertex_compsz, graph);

    for (std::string nodeinfo : vm["nodes"].as<std::vector<std::string>>())
    {
        Node &node = node_map[boost::add_vertex(graph)];
        std::stringstream ss(nodeinfo);
        std::string is_root, is_target, v, u;
        std::getline(std::getline(std::getline(std::getline(std::getline(ss, node.name, ','), is_root, ','), is_target, ','), v, ','), u, ',');
        node.is_root = std::stoi(is_root);
        node.is_target = std::stoi(is_target);
        try{node.ve = std::stod(v);}catch(...){node.ve = 0.0;}
        try{node.ue = std::stod(u);}catch(...){node.ue = 0.0;}
    }

    if (vm.count("shorts"))
        for (std::string info : vm["shorts"].as<std::vector<std::string>>())
        {
            std::stringstream ss(info);
            std::string name, tail, head, match, mismatch, v, u, vb, ub, T, min_score;
            std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(ss, name, ','), tail, ','), head, ','), match, ','), mismatch, ','), v, ','), u, ','), T, ','), min_score, ','), vb, ','), ub, ',');
            boost::graph_traits<graph_t>::edge_descriptor ed;
            bool inserted = false;
            for (boost::graph_traits<graph_t>::vertex_descriptor n = 0, t = boost::num_vertices(graph), h = boost::num_vertices(graph); n < boost::num_vertices(graph); ++n)
            {
                if (t == boost::num_vertices(graph) && node_map[n].name == tail)
                    t = n;
                if (h == boost::num_vertices(graph) && node_map[n].name == head)
                    h = n;
                if (t != boost::num_vertices(graph) && h != boost::num_vertices(graph))
                {
                    boost::tie(ed, inserted) = boost::add_edge(t, h, graph);
                    break;
                }
            }
            if (!inserted)
            {
                std::cerr << "short edge " << name << " cannot be added\n";
                exit(EXIT_FAILURE);
            }
            Edge &edge = edge_map[ed];
            edge.n = boost::num_edges(graph) - 1;
            edge.name = name;
            for (NUCTYPE a = 2; a < RankVec::sigma; ++a)
                for (NUCTYPE b = 2; b < RankVec::sigma; ++b)
                    if (a==b && a>2)
                        try{edge.gamma[a][b]=std::stod(match);}catch(...){edge.gamma[a][b]=1;}
                    else
                        try{edge.gamma[a][b]=std::stod(mismatch);}catch(...){edge.gamma[a][b]=-3;}
            try{edge.ve = std::stod(v);}catch(...){edge.ve = -5;}
            edge.vf = edge.ve;
            try{edge.ue = std::stod(u);}catch(...){edge.ue = -2;}
            edge.uf = edge.ue;
            try{edge.T = std::stod(T);}catch(...){edge.T = -10;}
            try{edge.min_score = std::stod(min_score);}catch(...){edge.min_score = 20;}
            try{edge.vfp = std::stod(vb);}catch(...){edge.vfp = 0;}
            edge.vfm = edge.vfp;
            try{edge.ufp = std::stod(ub);}catch(...){edge.ufp = 0;}
            edge.ufm = edge.ufp;

            SHORTSIZE ref_sz = file2short[edge.name].second;
            get_affine(edge.gf, ref_sz, 0, edge.uf, edge.vf);
            edge.gfm = edge.vfm + (ref_sz - 1) * edge.ufm;
            get_affine(edge.gfpT, ref_sz, edge.T, edge.ufp, edge.vfp);
        }

    if (vm.count("longs"))
        for (std::string info : vm["longs"].as<std::vector<std::string>>())
        {
            std::stringstream ss(info);
            std::string name, tail, head, match, mismatch, v, u, vb, ub, T, min_score;
            std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(ss, name, ','), tail, ','), head, ','), match, ','), mismatch, ','), v, ','), u, ','), T, ','), min_score, ',');
            boost::graph_traits<graph_t>::edge_descriptor ed;
            bool inserted = false;
            for (boost::graph_traits<graph_t>::vertex_descriptor n = 0, t = boost::num_vertices(graph), h = boost::num_vertices(graph); n < boost::num_vertices(graph); ++n)
            {
                if (t == boost::num_vertices(graph) && node_map[n].name == tail)
                    t = n;
                if (h == boost::num_vertices(graph) && node_map[n].name == head)
                    h = n;
                if (t != boost::num_vertices(graph) && h != boost::num_vertices(graph))
                {
                    boost::tie(ed, inserted) = boost::add_edge(t, h, graph);
                    break;
                }
            }
            if (!inserted)
            {
                std::cerr << "short edge " << name << " cannot be added\n";
                exit(EXIT_FAILURE);
            }
            Edge &edge = edge_map[ed];
            edge.n = boost::num_edges(graph) - 1;
            edge.name = name;

            for (NUCTYPE a = 2; a < RankVec::sigma; ++a)
                for (NUCTYPE b = 2; b < RankVec::sigma; ++b)
                    if (a==b && a>2)
                        try{edge.gamma[a][b]=std::stod(match);}catch(...){edge.gamma[a][b]=1;}
                    else
                        try{edge.gamma[a][b]=std::stod(mismatch);}catch(...){edge.gamma[a][b]=-3;}
            try{edge.ve = std::stod(v);}catch(...){edge.ve = -5;}
            edge.vf = edge.ve;
            try{edge.ue = std::stod(u);}catch(...){edge.ue = -2;}
            edge.uf = edge.ue;
            try{edge.T = std::stod(T);}catch(...){edge.T = -10;}
            try{edge.min_score = std::stod(min_score);}catch(...){edge.min_score = 20;}
        }
    
    std::vector<SIZETYPE> t_map(boost::num_vertices(graph));
    SIZETYPE comp_num = boost::strong_components(graph, comp_map, boost::discover_time_map(boost::make_iterator_property_map(t_map.begin(), boost::get(boost::vertex_index, graph))));
    std::vector<SIZETYPE> comp_time(comp_num, std::numeric_limits<int>::max());
    for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
        comp_time[comp_map[nd]] = std::min(comp_time[comp_map[nd]], t_map[nd]);
    std::vector<SIZETYPE> idx(comp_num);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&comp_time](SIZETYPE idx1, SIZETYPE idx2){return comp_time[idx1] < comp_time[idx2];});
    for (SIZETYPE i = 0; i < comp_num; ++i)
        comp_time[idx[i]] = i;
    for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
        comp_map[nd] = comp_time[comp_map[nd]];

    std::vector<SIZETYPE> scc_szs(comp_num, 0);
    for (boost::graph_traits<graph_t>::vertex_descriptor n = 0; n < boost::num_vertices(graph); ++n)
        ++scc_szs[comp_map[n]];
    for (boost::graph_traits<graph_t>::vertex_descriptor n = 0; n < boost::num_vertices(graph); ++n)
        compsz_map[n] = scc_szs[comp_map[n]];
    
    return comp_num;
}

int draw(graph_t &graph, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short)
{
    std::ofstream fout("graph.gv");
    fout << "digraph align_graph\n{\n";
    fout << "\tnode[shape=circle]\n";

    boost::property_map<graph_t, boost::vertex_Node_t>::type node_map = boost::get(boost::vertex_Node, graph);
    for (boost::graph_traits<graph_t>::vertex_descriptor n = 0; n < boost::num_vertices(graph); ++n)
    {
        Node &node = node_map[n];
        if (node.is_root)
            fout << '\t' << node.name << "[shape=square]\n";
        else if (node.is_target)
            fout << '\t' << node.name << "[shape=diamond]\n";
    }
    
    boost::property_map<graph_t, boost::edge_Edge_t>::type edge_map = boost::get(boost::edge_Edge, graph);
    boost::property_map<graph_t, boost::vertex_comp_t>::type comp_map = boost::get(boost::vertex_comp, graph);
    boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
    {
        Edge &edge = edge_map[*ei];
        boost::graph_traits<graph_t>::vertex_descriptor t = boost::source(*ei, graph), h = boost::target(*ei, graph); 
        Node &tail = node_map[t], &head = node_map[h];
        fout << '\t' << tail.name << "->" << head.name;
        if (file2short.count(edge.name))
        {
            if (comp_map[h] != comp_map[t])
                fout << "[color=black,label=\"";
            else
                fout << "[color=red,label=\"";
        }
        else
        {
            if (comp_map[h] != comp_map[t])
                fout << "[color=green,label=\"";
            else
                fout << "[color=blue,label=\"";
        }
        fout << edge.name << "\"]\n";
    }
    fout << "}\n";
    fout.close();
    return system("dot -Tpdf graph.gv > graph.gv.pdf");
}

#endif
