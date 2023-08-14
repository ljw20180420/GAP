#ifndef ALIGN_H
#define ALIGN_H

#include <climits>
#include <new>
#include <initializer_list>
#include <iomanip>
#include <cmath>
#include <set>
#include <typeinfo>
#include <bit>
#include "Graph.h"

struct Align
{
    graph_t graph;
    std::vector<SIZETYPE> comp_time;
    boost::property_map<graph_t, boost::vertex_Node_t>::type node_map;
    boost::property_map<graph_t, boost::edge_Edge_t>::type edge_map;
    boost::property_map<graph_t, boost::vertex_comp_t>::type comp_map;
    boost::property_map<graph_t, boost::vertex_compsz_t>::type compsz_map;
    std::vector<boost::graph_traits<graph_t>::edge_iterator> eis;

    std::mutex &mtx;
    std::ifstream &fin;
    std::vector<NUCTYPE> O;
    std::string Oname;
    std::ofstream fout;
    IDTYPE max_id = 0;
    QUERYSIZE Omax;
    std::unique_ptr<IDTYPE []> ids;
    std::unique_ptr<BITSTYPE []> bits;

    std::unique_ptr<SCORETYPE []> Ethres, croF;
    std::unique_ptr<QUERYSIZE []> crow;
    SIZETYPE crosize;

    SIZETYPE trn, tnn;

    std::unique_ptr<SCORETYPE []> vals;
    std::unique_ptr<SHORTSIZE []> ss;
    std::unique_ptr<DOTTYPE []> ns;

    SCORETYPE *pQval;

    std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short;
    std::map<std::string, RankVec> &file2rankvec;
    std::map<std::string, std::ifstream> file2long;
    std::map<std::string, std::ifstream> file2SA;

    Align(std::mutex &mtx_, std::ifstream &fin_, boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short_, std::map<std::string, RankVec> &file2rankvec_, std::string mgfile, QUERYSIZE Omax_)
        : mtx(mtx_), fin(fin_), file2short(file2short_), file2rankvec(file2rankvec_), fout(mgfile, std::ifstream::binary), Omax(Omax_)
    {
        comp_time = construct_graph(graph, vm, file2short);
        node_map = boost::get(boost::vertex_Node, graph);
        edge_map = boost::get(boost::edge_Edge, graph);
        comp_map = boost::get(boost::vertex_comp, graph);
        compsz_map = boost::get(boost::vertex_compsz, graph);
        eis.resize(boost::num_edges(graph));
        boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
            eis[edge_map[*ei].n] = ei;

        if (vm.count("longs"))
            for (std::string file : vm["longs"].as<std::vector<std::string>>())
            {
                file = file.substr(0, file.find(','));
                if (file2long.count(file))
                    continue;
                file2long[file].open(file);
                file2SA[file].open(file+".sa");
            }

        trn = 0;
        tnn = 1 + boost::num_vertices(graph); // Q and Abars
        for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
        {
            trn += compsz_map[nd] + 1;
            tnn += (Omax + 1) * (compsz_map[nd] + 1);
        }

        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (comp_map[boost::target(*ei, graph)] == comp_map[boost::source(*ei, graph)] || !file2short.count(edge.name))
                continue;
            SHORTSIZE ref_sz = file2short[edge.name].second;
            trn += 3 * (ref_sz + 1);
            tnn += (Omax + 1) * 3 * (ref_sz + 1);
        }
        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (comp_map[boost::target(*ei, graph)] != comp_map[boost::source(*ei, graph)] || !file2short.count(edge.name))
                continue;
            SIZETYPE tail_scc_sz = compsz_map[boost::source(*ei, graph)];
            SHORTSIZE ref_sz = file2short[edge.name].second;
            trn += 1 + 4 * (ref_sz + 1) + 2 * (tail_scc_sz - 1);
            tnn += (Omax + 1) * (1 + 4 * (ref_sz + 1) + 2 * (tail_scc_sz - 1));
        }
        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (comp_map[boost::target(*ei, graph)] != comp_map[boost::source(*ei, graph)] || file2short.count(edge.name))
                continue;
            SIZETYPE tail_scc_sz = compsz_map[boost::source(*ei, graph)];
            trn += tail_scc_sz - 1;
            tnn += (Omax + 1) * (tail_scc_sz - 1);
        }
        vals.reset(new SCORETYPE[tnn]);
        ids.reset(new IDTYPE[tnn]);
        ss.reset(new SHORTSIZE[trn]);
        ns.reset(new DOTTYPE[trn]);
        bits.reset(new BITSTYPE[tnn]);

        SCORETYPE *fpval = vals.get(), *rpval = fpval + tnn;
        pQval = --rpval;
        SHORTSIZE *fps = ss.get();
        DOTTYPE *fpn = ns.get();
        for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
        {
            Node &node = node_map[nd];
            node.pAbarval = --rpval;
            *(fps++) = 0;
            *(fpn++) = Dot::NODEIDX2DOTTYPEB(nd);
            for (SIZETYPE j = 0; j < compsz_map[nd]; ++j, ++fps, ++fpn)
            {
                *fps = j;
                *fpn = Dot::NODEIDX2DOTTYPEA(nd);
            }
            node.Bvals = fpval;
            fpval += Omax + 1;
            node.Avals = new SCORETYPE *[compsz_map[nd]];
            for (SIZETYPE s = 0; s < compsz_map[nd]; ++s, fpval += Omax + 1)
                node.Avals[s] = fpval;
            node.AdeltaDot = new std::deque<SCORETYPE *>[Omax + 1];
            node.AdeltaGlobal = new std::deque<Node::GlobalSuffix>[Omax + 1];
        }

        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (comp_map[boost::target(*ei, graph)] == comp_map[boost::source(*ei, graph)] || !file2short.count(edge.name))
                continue;
            SHORTSIZE ref_sz = file2short[edge.name].second;
            for (SIZETYPE i = 0; i < 3; ++i)
                for (SIZETYPE j = 0; j <= ref_sz; ++j, ++fps, ++fpn)
                {
                    *fps = j;
                    *fpn = edge.n;
                }
            SIZETYPE si = 0;
            for (SCORETYPE ***pYvalss : {&edge.Evals, &edge.Fvals, &edge.Gvals})
            {
                *pYvalss = new SCORETYPE *[ref_sz + 1];
                for (SIZETYPE s = 0; s <= ref_sz; ++s, fpval += Omax + 1)
                    (*pYvalss)[s] = fpval;
            }
        }

        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (comp_map[boost::target(*ei, graph)] != comp_map[boost::source(*ei, graph)] || !file2short.count(edge.name))
                continue;
            SHORTSIZE ref_sz = file2short[edge.name].second;
            *(fps++) = 0;
            *(fpn++) = edge.n;
            SIZETYPE tail_scc_sz = compsz_map[boost::source(*ei, graph)];
            SIZETYPE Ss[6] = {ref_sz+1, ref_sz+1, ref_sz+1, ref_sz+1, tail_scc_sz-1, tail_scc_sz-1};
            for (SIZETYPE i = 0; i < 6; ++i)
                for (SIZETYPE j = 0; j < Ss[i]; ++j, ++fps, ++fpn)
                {
                    *fps = j;
                    *fpn = edge.n;
                }
            edge.Dvals = fpval;
            fpval += Omax + 1;
            SIZETYPE si = 0;
            for (SCORETYPE ***pYvalss : {&edge.Evals, &edge.F0vals, &edge.G0vals, &edge.Gvals, &edge.D0vals, &edge.DXvals})
            {
                *pYvalss = new SCORETYPE *[Ss[si]];
                for (SIZETYPE s = 0; s < Ss[si]; ++s, fpval += Omax + 1)
                    (*pYvalss)[s] = fpval;
                ++si;
            }
        }

        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (comp_map[boost::target(*ei, graph)] != comp_map[boost::source(*ei, graph)] || file2short.count(edge.name))
                continue;
            SIZETYPE tail_scc_sz = compsz_map[boost::source(*ei, graph)];
            for (SIZETYPE j = 0; j < tail_scc_sz-1; ++j, ++fps, ++fpn)
            {
                *fps = j;
                *fpn = edge.n;
            }
            edge.D0vals = new SCORETYPE *[tail_scc_sz-1];
            for (SIZETYPE s = 0; s < tail_scc_sz-1; ++s, fpval += Omax + 1)
                edge.D0vals[s] = fpval;
        }

        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (!file2short.count(edge.name))
                edge.tailAvals = new SCORETYPE[Omax + 1];
        }

        Ethres.reset(new SCORETYPE[Omax + 1]);
        crosize = std::ceil((Omax + 1) * (Omax + 1) * 1.4); // the initial crosize assume matchscore = 1 and u = -2, so the aligning shift is upper-bounded by 1/3, thereby 1.4 > 1 + 1/3 should be safe
        croF.reset(new SCORETYPE[2 * crosize]);
        crow.reset(new QUERYSIZE[crosize]);
    }

    void run()
    {
        std::string Otmp;
        while (true)
        {
            mtx.lock();
            if (!std::getline(std::getline(fin, Oname), Otmp))
            {
                mtx.unlock();
                break;
            }
            mtx.unlock();
            O.clear();
            for (SIZETYPE i = 0; i < Otmp.size(); ++i)
                O.push_back(base2NUCTYPE[Otmp[i]]);
            Mix();
            GetMinimalGraph();
        }
        if (max_id>0)
            fout.write((char *)&max_id, sizeof(max_id));
        fout.close(); // close is must because even after run() return, fout as a member of align, is not destructed, thereby its buffer is not output
    }

    void Mix()
    {
        for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
        {
            Node &node = node_map[nd];
            for (SIZETYPE w = 0; w <= O.size(); ++w)
                node.Avals[0][w] = -inf;
            if (node.is_root)
            {
                *node.pAbarval = 0;
                node.Avals[0][0] = 0;
                node.AdeltaDot[0].emplace_back(node.pAbarval);
            }
        }

        boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
        SIZETYPE scc_sz;
        for (SIZETYPE ct : comp_time)
        {
            for (SIZETYPE w = 0; w <= O.size(); ++w)
            {
                for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
                {
                    if (comp_map[nd] != ct)
                        continue;

                    Node &node = node_map[nd];
                    SIZETYPE M = &node.Bvals[w] - vals.get();
                    bits[M] = 0;
                    if (w == 0)
                        node.Bvals[w] = -inf;
                    else
                    {
                        node.Bvals[w] = std::max(node.Bvals[w - 1] + node.ue, node.Avals[compsz_map[nd] - 1][w - 1] + node.ve);
                        if (node.Bvals[w] == node.Bvals[w - 1] + node.ue)
                            bits[M] += 1;
                        if (node.Bvals[w] == node.Avals[compsz_map[nd] - 1][w - 1] + node.ve)
                            bits[M] += 2;                    
                    }
                    if (node.Bvals[w] >= node.Avals[0][w])
                    {
                        if (node.Bvals[w] > node.Avals[0][w])
                        {
                            node.Avals[0][w] = node.Bvals[w];
                            node.AdeltaDot[w].clear();
                            node.AdeltaGlobal[w].clear();
                        }
                        node.AdeltaDot[w].emplace_back(&node.Bvals[w]);
                    }    

                    scc_sz = compsz_map[nd];
                }
                for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
                    if (comp_map[boost::target(*ei, graph)] == ct && comp_map[boost::source(*ei, graph)] == ct)
                    {
                        if (file2short.count(edge_map[*ei].name))
                            CircuitIteration(w, ei);
                        else
                            CircuitIterationGlobal(w, ei);
                    }
                for (SIZETYPE l = 1; l < scc_sz; ++l)
                {
                    for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
                    {
                        Edge &edge = edge_map[*ei];
                        if (comp_map[boost::target(*ei, graph)] != ct || comp_map[boost::source(*ei, graph)] != ct || !file2long.count(edge.name))
                            continue;
                        edge.D0vals[l - 1][w] = node_map[boost::source(*ei, graph)].Avals[l - 1][w] + edge.T;
                    }
                    for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
                    {
                        Edge &edge = edge_map[*ei];
                        if (comp_map[boost::target(*ei, graph)] != ct || comp_map[boost::source(*ei, graph)] != ct || !file2short.count(edge.name))
                            continue;
                        Node &tail = node_map[boost::source(*ei, graph)];
                        SHORTSIZE ref_sz = file2short[edge.name].second;
                        edge.D0vals[l - 1][w] = tail.Avals[l - 1][w] + edge.T;
                        SIZETYPE M = &(edge.DXvals[l - 1][w]) - vals.get();
                        bits[M] = 0;
                        edge.DXvals[l - 1][w] = std::max(tail.Avals[l - 1][w] + edge.gfpT[ref_sz], edge.D0vals[l - 1][w] + edge.gf[ref_sz]);
                        if (edge.DXvals[l - 1][w] == tail.Avals[l - 1][w] + edge.gfpT[ref_sz])
                            bits[M] += 1;
                        if (edge.DXvals[l - 1][w] == edge.D0vals[l - 1][w] + edge.gf[ref_sz])
                            bits[M] += 2;
                    }
                    for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
                    {
                        if (comp_map[nd] != ct)
                            continue;
                        
                        Node &node = node_map[nd];
                        SIZETYPE M = &node.Avals[l][w] - vals.get();
                        node.Avals[l][w] = node.Avals[0][w];
                        bits[M] = 1;
                        SIZETYPE ib = 1;
                        boost::graph_traits<graph_t>::in_edge_iterator inei, inei_end;
                        for (boost::tie(inei, inei_end) = boost::in_edges(nd, graph); inei != inei_end; ++inei)
                        {
                            Edge &edge = edge_map[*inei];
                            if (comp_map[boost::source(*inei, graph)] != ct)
                                continue;
                            ib *= 2;
                            if (ib > std::numeric_limits<BITSTYPE>::max())
                            {
                                std::cerr << "the upper limits of BITSTYPE reached by node" << nd << ", A[" << l << "][" << w << "]\n";
                                break;
                            }
                            SCORETYPE Dscore;
                            if (file2short.count(edge.name))
                                Dscore = std::max(edge.D0vals[l - 1][w] + edge.gfm, edge.DXvals[l - 1][w]);
                            else
                                Dscore = edge.D0vals[l - 1][w];
                            if (Dscore >= node.Avals[l][w])
                            {
                                if (Dscore > node.Avals[l][w])
                                {
                                    node.Avals[l][w] = Dscore;
                                    bits[M] = 0;
                                }
                                bits[M] += ib;
                            }
                        }
                    }
                }
            }
            for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
                if (comp_map[boost::source(*ei, graph)] == ct && comp_map[boost::target(*ei, graph)] != ct && file2short.count(edge_map[*ei].name))
                    CrossIteration(ei);
            for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
                if (comp_map[boost::source(*ei, graph)] == ct && comp_map[boost::target(*ei, graph)] != ct && file2long.count(edge_map[*ei].name))
                    CrossIterationGlobal(ei);
        }

        SIZETYPE M = pQval - vals.get();
        bits[M] = 0;
        *pQval = -inf;
        SIZETYPE ib = 1;
        for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
        {
            Node &node = node_map[nd];
            if (!node.is_target)
                continue;
            if (ib > std::numeric_limits<BITSTYPE>::max())
            {
                std::cerr << "the upper limits of BITSTYPE reached by Q (the graph has too many targets)\n";
                break;
            }
            SCORETYPE score = node.Avals[compsz_map[nd] - 1][O.size()];
            if (score >= *pQval)
            {
                if (score > *pQval)
                {
                    *pQval = score;
                    bits[M] = 0;
                }
                bits[M] += ib;
            }
            ib *= 2;
        }
    }

    void CrossIteration(boost::graph_traits<graph_t>::edge_iterator ei)
    {
        Edge &edge = edge_map[*ei];
        NUCTYPE *ref = file2short[edge.name].first.get();
        SHORTSIZE ref_sz = file2short[edge.name].second;
        SCORETYPE *Avals = node_map[boost::source(*ei, graph)].Avals[compsz_map[boost::source(*ei, graph)] - 1];
        SCORETYPE **Evals = edge.Evals, **Fvals = edge.Fvals, **Gvals = edge.Gvals;
        for (SIZETYPE s = 0; s <= ref_sz; ++s)
            for (SIZETYPE w = 0; w <= O.size(); ++w)
            {
                SIZETYPE M = &Evals[s][w] - vals.get();
                bits[M] = 0;
                if (w == 0)
                    vals[M] = -inf;
                else
                {
                    vals[M] = std::max(Evals[s][w - 1] + edge.ue, Gvals[s][w - 1] + edge.ve);
                    if (vals[M] == Evals[s][w - 1] + edge.ue)
                        bits[M] += 1;
                    if (vals[M] == Gvals[s][w - 1] + edge.ve)
                        bits[M] += 2;
                }

                M = &Fvals[s][w] - vals.get();
                bits[M] = 0;
                if (s == 0)
                    vals[M] = -inf;
                else
                {
                    vals[M] = std::max(Fvals[s - 1][w] + edge.uf, Gvals[s - 1][w] + edge.vf);
                    if (vals[M] == Fvals[s - 1][w] + edge.uf)
                        bits[M] += 1;
                    if (vals[M] == Gvals[s - 1][w] + edge.vf)
                        bits[M] += 2;
                }
                
                SCORETYPE scores[4] = {Evals[s][w], Fvals[s][w], Avals[w] + edge.gfpT[s], -inf};
                if (s > 0 && w > 0)
                    scores[3] = Gvals[s - 1][w - 1] + edge.gamma[ref[s - 1]][O[w - 1]];
                M = &Gvals[s][w] - vals.get();
                vals[M] = scores[0];
                bits[M] = 1;
                for (SIZETYPE i = 1, bi = 2; i < 4; ++i, bi *= 2)
                    if (scores[i] >= vals[M])
                    {
                        if (scores[i] > vals[M])
                        {
                            vals[M] = scores[i];
                            bits[M] = 0;
                        }
                        bits[M] += bi;
                    }

                Node &head = node_map[boost::target(*ei, graph)];
                if (Gvals[s][w] >= head.Avals[0][w])
                {
                    if (Gvals[s][w] > head.Avals[0][w])
                    {
                        head.Avals[0][w] = Gvals[s][w];
                        head.AdeltaDot[w].clear();
                        head.AdeltaGlobal[w].clear();
                    }
                    head.AdeltaDot[w].emplace_back(&Gvals[s][w]);
                }
            }
    }

    struct SimNode
    {
        NUCTYPE c;
        SIZETYPE sr1, sr2;
        SIZETYPE shiftFG;

        SimNode(NUCTYPE c_, SIZETYPE sr1_, SIZETYPE sr2_, SIZETYPE shiftFG_)
        : c(c_), sr1(sr1_), sr2(sr2_), shiftFG(shiftFG_)
        {}
    };

    void CrossIterationGlobal(boost::graph_traits<graph_t>::edge_iterator ei)
    {
        Edge &edge = edge_map[*ei];
        Node &head = node_map[boost::target(*ei, graph)], &tail = node_map[boost::source(*ei, graph)];
        SCORETYPE *Avals = tail.Avals[compsz_map[boost::source(*ei, graph)] - 1];
        RankVec &rankvec = file2rankvec[edge.name];
        SCORETYPE *croG = croF.get() + crosize;

        std::deque<SimNode> simnodes = {SimNode{0, 0, rankvec.bwt_sz, 0}}; // simnodes[0] has no base, we simply set it to 0, which does not mean that it has base #(0)
        SCORETYPE EO, EN;
        for (SIZETYPE w = 0; w <= O.size(); ++w)
        {
            if (w > 0)
                EN = std::max(EO + edge.ue, croG[w - 1] + edge.ve);
            else
                EN = -inf;
            SCORETYPE tgw = (O.size() > w ? (O.size() - w - 1) * tail.ue + tail.ve : 0);
            crow[w] = w;
            croF[w] = -inf;
            croG[w] = std::max(EN, Avals[w] + edge.T);
            Ethres[w] = std::max(EN, croG[w] + std::min({SCORETYPE(0), tail.ve - edge.ue, SCORETYPE(tgw - (O.size() - w) * edge.ue)}));

            if (croG[w] >= head.Avals[0][w])
            {
                if (croG[w] > head.Avals[0][w])
                {
                    head.Avals[0][w] = croG[w];
                    head.AdeltaDot[w].clear();
                    head.AdeltaGlobal[w].clear();
                }
                head.AdeltaGlobal[w].emplace_back(ei, simnodes[0].sr1, 0);
            }

            EO = EN;
        }
        SIZETYPE shiftFG = O.size() + 1;


        while (true)
        {
            while (true)
            {
                NUCTYPE c;
                if (simnodes.back().shiftFG < shiftFG)
                    c = 1;
                else
                {
                    do
                    {
                        if (simnodes.size() == 1)
                            return;
                        shiftFG = simnodes.back().shiftFG;
                        c = simnodes.back().c;
                        simnodes.pop_back();
                    } while(c == 6);
                }
                SIZETYPE sr1, sr2;
                do
                {
                    ++c;
                    sr1 = rankvec.C[c] + rankvec.rank(c, simnodes.back().sr1);
                    sr2 = rankvec.C[c] + rankvec.rank(c, simnodes.back().sr2);
                } while (sr1 >= sr2 && c < 6);
                if (sr1 < sr2)
                {
                    simnodes.emplace_back(c, sr1, sr2, shiftFG);
                    break;
                }
                else 
                    shiftFG = simnodes.back().shiftFG;
            }

            for (SIZETYPE idx = simnodes[simnodes.size() - 2].shiftFG, w; idx < simnodes.back().shiftFG; ++idx)
            {
                if (idx >= crosize)
                {
                    std::cerr << "the default memory is too low for CrossIterationGlobal\n";
                    exit(EXIT_FAILURE);
                }
                if (idx == simnodes[simnodes.size() - 2].shiftFG || w < crow[idx] - 1) // if w exit loop not by break, then w = crow[idx]
                    EO = -inf;
                else
                    EO = EN;
                QUERYSIZE W;
                if (idx < simnodes.back().shiftFG - 1)
                    W = crow[idx + 1];
                else
                    W = O.size() + 1;
                for (w = crow[idx]; w < W; ++w)
                {
                    if (w == 0)
                        EN = -inf;
                    else
                        EN = std::max(EO + edge.ue, (simnodes.back().shiftFG < shiftFG && crow[shiftFG - 1] == w - 1) ? croG[shiftFG - 1] + edge.ve : -inf);
                    if (EN < Ethres[w])
                        EN = -inf;
                    SCORETYPE F_val;
                    if (w == crow[idx])
                        F_val = std::max(croF[idx] + edge.uf, croG[idx] + edge.vf);
                    else
                        F_val = -inf;
                    SCORETYPE G_val;
                    if (w == crow[idx] && idx > simnodes[simnodes.size() - 2].shiftFG && crow[idx - 1] == w - 1)
                        G_val = croG[idx - 1] + edge.gamma[simnodes.back().c][O[w - 1]];
                    else if (w == crow[idx] + 1)
                        G_val = croG[idx] + edge.gamma[simnodes.back().c][O[w - 1]];
                    else
                        G_val = -inf;
                    G_val = std::max({G_val, EN, F_val});
                    if (G_val >= croG[w])
                    {
                        crow[shiftFG] = w;
                        croF[shiftFG] = F_val;
                        croG[shiftFG] = G_val;
                        ++shiftFG;

                        if (G_val >= head.Avals[0][w])
                        {
                            if (G_val > head.Avals[0][w])
                            {
                                head.Avals[0][w] = G_val;
                                head.AdeltaDot[w].clear();
                                head.AdeltaGlobal[w].clear();
                            }
                            head.AdeltaGlobal[w].emplace_back(ei, simnodes.back().sr1, simnodes.size() - 1);
                        }
                    }
                    else if (w > crow[idx] + 1 && EN == -inf)
                        break;
                    EO = EN;
                }
            }
        }
    }

    void CircuitIteration(QUERYSIZE w, boost::graph_traits<graph_t>::edge_iterator ei)
    {
        Edge &edge = edge_map[*ei];
        Node &head = node_map[boost::target(*ei, graph)], &tail = node_map[boost::source(*ei, graph)];
        NUCTYPE *ref = file2short[edge.name].first.get();
        SHORTSIZE ref_sz = file2short[edge.name].second;
        SCORETYPE *Avals = tail.Avals[compsz_map[boost::source(*ei, graph)] - 1];
        SCORETYPE *Dvals = edge.Dvals;
        SCORETYPE **Evals = edge.Evals, **F0vals = edge.F0vals, **G0vals = edge.G0vals, **Gvals = edge.Gvals;

        if (w > 0)
            Dvals[w - 1] = Avals[w - 1] + edge.T;

        for (SIZETYPE s = 0; s <= ref_sz; ++s)
        {
            if (w > 0)
            {
                SCORETYPE scores[3] = {G0vals[s][w - 1], Dvals[w - 1] + edge.gf[s], -inf};
                if (s > 0)
                    scores[2] = Avals[w - 1] + edge.gfpT[s];
                SIZETYPE M = &(Gvals[s][w - 1]) - vals.get();
                vals[M] = scores[0];
                bits[M] = 1;
                for (SIZETYPE i = 1, bi = 2; i < 3; ++i, bi *= 2)
                {
                    if (scores[i] >= vals[M])
                    {
                        if (scores[i] > vals[M])
                        {
                            vals[M] = scores[i];
                            bits[M] = 0;
                        }
                        bits[M] += bi;
                    }
                }
            }

            SIZETYPE M = &Evals[s][w] - vals.get();
            bits[M] = 0;
            if (w == 0)
                vals[M] = -inf;
            else
            {
                vals[M] = std::max(Evals[s][w - 1] + edge.ue, Gvals[s][w - 1] + edge.ve);
                if (vals[M] == Evals[s][w - 1] + edge.ue)
                    bits[M] += 1;
                if (vals[M] == Gvals[s][w - 1] + edge.ve)
                    bits[M] += 2;
            }
                
            M = &F0vals[s][w] - vals.get();
            bits[M] = 0;
            if (s == 0)
                vals[M] = -inf;
            else
            {
                vals[M] = std::max(F0vals[s - 1][w] + edge.uf, G0vals[s - 1][w] + edge.vf);
                if (vals[M] == F0vals[s - 1][w] + edge.uf)
                    bits[M] += 1;
                if (vals[M] == G0vals[s - 1][w] + edge.vf)
                    bits[M] += 2;
            }

            SCORETYPE scores[3] = {Evals[s][w], F0vals[s][w], -inf};
            if (s > 0 && w > 0)
                scores[2] = Gvals[s - 1][w - 1] + edge.gamma[ref[s - 1]][O[w - 1]];
            M = &G0vals[s][w] - vals.get();
            vals[M] = scores[0];
            bits[M] = 1;
            for (SIZETYPE i = 1, bi = 2; i < 3; ++i, bi *= 2)
            {
                if (scores[i] >= vals[M])
                {
                    if (scores[i] > vals[M])
                    {
                        vals[M] = scores[i];
                        bits[M] = 0;
                    }
                    bits[M] += bi;
                }
            }

            if (G0vals[s][w] >= head.Avals[0][w])
            {
                if (G0vals[s][w] > head.Avals[0][w])
                {
                    head.Avals[0][w] = G0vals[s][w];
                    head.AdeltaDot[w].clear();
                    head.AdeltaGlobal[w].clear();
                }
                head.AdeltaDot[w].emplace_back(&G0vals[s][w]);
            }
        }
    }

    void CircuitIterationGlobal(QUERYSIZE w, boost::graph_traits<graph_t>::edge_iterator ei)
    {
        if (w == 0)
            return;

        SCORETYPE tgws;
        SNC *jump;
        RankVec *prankvec = &file2rankvec[edge_map[*ei].name];

        boost::graph_traits<graph_t>::vertex_descriptor hd = boost::target(*ei, graph), td = boost::source(*ei, graph);
        Edge &edge = edge_map[*ei];
        Node &head = node_map[hd], &tail = node_map[td];
        if (w == 1)
        {
            edge.sncs.emplace_back(-inf, -inf, -inf, -inf, 0, 5, 0, 0, prankvec->bwt_sz, 1, (SNC *)NULL, (SNC *)NULL); // the root SNC has no base, we simply set it to 0, which does not mean that it has base #(0)
            for (NUCTYPE c = 2; c <= 6; ++c)
            {
                SIZETYPE sr1 = prankvec->C[c] + prankvec->rank(c, edge.sncs.front().sr1);
                SIZETYPE sr2 = prankvec->C[c] + prankvec->rank(c, edge.sncs.front().sr2);
                if (sr1 < sr2)
                {
                    edge.sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, 1, sr1, sr2, 0, &(edge.sncs.front()), edge.sncs.front().jump);
                    edge.sncs.front().jump = &(edge.sncs.back());
                }
            }
        }
        else
        {
            for (SNC *prejump = &(edge.sncs.front()), *jump = prejump->jump; jump; jump = prejump->jump)
            {
                jump->hatG = jump->G0;
                if (jump->itp != &(edge.sncs.front()) && jump->itp->hatG == -inf && jump->hatG == -inf && jump->E == -inf)
                    prejump->jump = jump->jump;
                else
                    prejump = jump;
            }
        }

        edge.sncs.front().hatG = std::max(edge.sncs.front().G0, tail.Avals[compsz_map[td] - 1][w - 1] + edge.T);
        edge.sncs.front().E = std::max(edge.sncs.front().hatG + edge.ve, edge.sncs.front().E + edge.ue);
        tgws = (O.size() > w ? (O.size() - w - 1) * tail.ue + tail.ve : 0);
        edge.sncs.front().F0 = -inf;
        edge.sncs.front().G0 = edge.sncs.front().E;

        if (edge.sncs.front().G0 >= head.Avals[0][w])
        {
            if (edge.sncs.front().G0 > head.Avals[0][w])
            {
                head.Avals[0][w] = edge.sncs.front().G0;
                head.AdeltaDot[w].clear();
                head.AdeltaGlobal[w].clear();
            }
            head.AdeltaGlobal[w].emplace_back(ei, edge.sncs.front().sr1, edge.sncs.front().lambda);
        }

        jump = edge.sncs.front().jump;

        while (true)
        {
            SCORETYPE Gthres = std::max(edge.sncs.front().E, tail.Avals[0][w] + edge.T);
            SCORETYPE Ethres = std::max(edge.sncs.front().E, Gthres + std::min(SCORETYPE(0), std::min(tail.ve - edge.ue, SCORETYPE(tgws - (O.size() - w) * edge.ue))));

            if (jump)
            {
                jump->E = std::max(jump->hatG + edge.ve, jump->E + edge.ue);
                if (jump->E < Ethres)
                    jump->E = -inf;
                jump->F0 = std::max(jump->itp->G0 + edge.vf, jump->itp->F0 + edge.uf);
                jump->G0 = std::max(std::max(jump->E, jump->F0), jump->itp->hatG + edge.gamma[jump->c][O[w - 1]]);
                if (jump->G0 < Gthres)
                {
                    jump->F0 = -inf;
                    jump->G0 = -inf;
                }
                else
                {
                    if (jump->G0 >= head.Avals[0][w])
                    {
                        if (jump->G0 > head.Avals[0][w])
                        {
                            head.Avals[0][w] = jump->G0;
                            head.AdeltaDot[w].clear();
                            head.AdeltaGlobal[w].clear();
                        }
                        head.AdeltaGlobal[w].emplace_back(ei, jump->sr1, jump->lambda);
                    }
                }
                    
                if (jump->G0 != -inf && jump->hatG == -inf)
                {
                    if (jump->cid == 0)
                    {
                        jump->cid = edge.sncs.size();
                        for (NUCTYPE c = 2; c < RankVec::sigma; ++c)
                        {
                            SIZETYPE sr1 = prankvec->C[c] + prankvec->rank(c, jump->sr1);
                            SIZETYPE sr2 = prankvec->C[c] + prankvec->rank(c, jump->sr2);
                            if (sr1 < sr2)
                            {
                                edge.sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, jump->lambda + 1, sr1, sr2, 0, jump, (SNC *)NULL);
                                ++jump->idl;
                            }
                        }
                    }
                    for (NUCTYPE idi = 0; idi < jump->idl; ++idi)
                        if (edge.sncs[jump->cid + idi].E == -inf && edge.sncs[jump->cid + idi].G0 == -inf)
                        {
                            edge.sncs[jump->cid + idi].jump = jump->jump;
                            jump->jump = &(edge.sncs[jump->cid + idi]);
                        }
                }
                jump = jump->jump;
            }
            else
            {
                if (w == O.size())
                    edge.sncs.clear();
                break;
            }
        }
    }

    struct SWN
    {
        SIZETYPE s;
        QUERYSIZE w;
        DOTTYPE n;
    };

    SWN get_swn(SIZETYPE idx)
    {
        SWN swn;
        swn.s = idx / (Omax + 1);
        if (swn.s < trn)
        {
            swn.n = ns[swn.s];
            swn.s = ss[swn.s];
            swn.w = idx % (Omax + 1);
        }
        else
        {
            swn.s = 0;
            swn.w = 0;
            if (idx == tnn - 1)
                swn.n = Dot::DotQ;
            else
                swn.n = Dot::DotAbar;
        }
        return swn;
    }

    enum VALTYPE{NA, NB, SOE, SOF, SOG, SID, SIE, SIF0, SIG0, SIG, SID0, SIDX, LID0, ABAR, Q};

    Align::VALTYPE get_type(SIZETYPE M, SWN &swn)
    {
        if (swn.n == Dot::DotQ)
            return VALTYPE::Q;
        if (swn.n == Dot::DotAbar)
            return VALTYPE::ABAR;
        if (swn.n>=0)
        {
            Edge &edge = edge_map[*eis[swn.n]];
            if (file2long.count(edge.name))
                return VALTYPE::LID0;
            if (comp_map[boost::source(*eis[swn.n], graph)] == comp_map[boost::target(*eis[swn.n], graph)])
            {
                if (M < edge.Evals[0] - vals.get())
                    return VALTYPE::SID;
                if (M < edge.F0vals[0] - vals.get())
                    return VALTYPE::SIE;
                if (M < edge.G0vals[0] - vals.get())
                    return VALTYPE::SIF0;
                if (M < edge.Gvals[0] - vals.get())
                    return VALTYPE::SIG0;
                if (M < edge.D0vals[0]-vals.get())
                    return VALTYPE::SIG;
                if (M < edge.DXvals[0]-vals.get())
                    return VALTYPE::SID0;
                return VALTYPE::SIDX;
            }
            else
            {
                if (M < edge.Fvals[0] - vals.get())
                    return VALTYPE::SOE;
                if (M < edge.Gvals[0] - vals.get())
                    return VALTYPE::SOF;
                return VALTYPE::SOG;
            }
        }
        if (swn.n % 2)
            return VALTYPE::NB;
        else
            return VALTYPE::NA;
    }

    void outputdot(SCORETYPE Mval, SIZETYPE M, SWN &swn, Align::VALTYPE type)
    {
        fout.write((char *)&swn.n, sizeof(Dot::n));
        fout.write((char *)&swn.s, sizeof(Dot::s));
        fout.write((char *)&swn.w, sizeof(Dot::w));
        fout.write((char *)&Mval, sizeof(Dot::val));
        SOURCESIZE s_sz;
        switch (type)
        {
        case VALTYPE::NB: case VALTYPE::SOE: case VALTYPE::SOF: case VALTYPE::SOG: case VALTYPE::SIDX: case VALTYPE::SIE: case VALTYPE::SIF0: case VALTYPE::SIG0: case VALTYPE::SIG: case VALTYPE::Q:
            s_sz = std::popcount(bits[M]);
            break;
        case VALTYPE::ABAR:
            s_sz = 0;
            break;
        case VALTYPE::LID0: case VALTYPE::SID0: case VALTYPE::SID:
            s_sz = 1;
            break;
        case VALTYPE::NA:
            if (swn.s == 0)
            {
                Node &node = node_map[Dot::DOTTYPE2NODEIDX(swn.n)];
                s_sz = node.AdeltaGlobal[swn.w].size() + node.AdeltaDot[swn.w].size();
            }
            else
                s_sz = std::popcount(bits[M]);
            break;
        }
        fout.write((char *)&s_sz, sizeof(Dot::s_sz));
        fout.write((char *)&s_sz, sizeof(Dot::lambda)); // strange
    }

    void outputdot(SCORETYPE Mval, Dot *M)
    {
        fout.write((char *)&M->n, sizeof(Dot::n));
        fout.write((char *)&M->s, sizeof(Dot::s));
        fout.write((char *)&M->w, sizeof(Dot::w));
        fout.write((char *)&Mval, sizeof(Dot::val));
        fout.write((char *)&M->s_sz, sizeof(Dot::s_sz));
        fout.write((char *)&M->lambda, sizeof(Dot::lambda));
    }

    void arrangesource(std::queue<SCORETYPE> &valuequeue, SCORETYPE *source, std::queue<SIZETYPE> &localqueue, IDTYPE &id)
    {
        SIZETYPE idx = source - vals.get();
        if (*source != -inf)
        {
            ids[idx] = ++id;
            localqueue.push(idx);
            valuequeue.push(*source);
            *source = -inf;
        }
        fout.write((char *)&ids[idx], sizeof(Dot::id));
    }

    void arrangesource(std::queue<SCORETYPE> &valuequeue, Dot *source, std::queue<Dot *> &globalqueue, IDTYPE &id)
    {
        if (source->val != -inf)
        {
            source->id = ++id;
            globalqueue.push(source);
            valuequeue.push(source->val);
            source->val = -inf;
        }
        fout.write((char *)&source->id, sizeof(Dot::id));
    }

    

    void GetMinimalGraph()
    {
        boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei)
        {
            Edge &edge = edge_map[*ei];
            if (file2long.count(edge.name))
            {
                edge.tracktree.initialize(edge.n, O.size());
                boost::graph_traits<graph_t>::vertex_descriptor td = boost::source(*ei, graph);
                SCORETYPE *Avals = node_map[td].Avals[compsz_map[td] - 1];
                for (SIZETYPE i=0; i<=O.size(); ++i)
                    edge.tailAvals[i] = Avals[i];
            }
        }

        SIZETYPE Onsz = Oname.size();
        fout.write((char *)&Onsz, sizeof(Onsz));
        fout.write(Oname.data(), Onsz);
        std::queue<Dot *> globalqueue;
        std::queue<SIZETYPE> localqueue;
        std::queue<SCORETYPE> valuequeue;
        IDTYPE id = 0;
        ids[pQval-vals.get()] = id;
        localqueue.push(pQval - vals.get());
        valuequeue.push(*pQval);
        *pQval = -inf;
        std::deque<Dot *> dot_sources;
        while (!globalqueue.empty() || !localqueue.empty())
        {
            if (globalqueue.empty() || !localqueue.empty() && ids[localqueue.front()] < globalqueue.front()->id)
            {
                SIZETYPE M = localqueue.front();
                localqueue.pop();
                SCORETYPE Mval = valuequeue.front();
                valuequeue.pop();
                SWN swn = get_swn(M);
                Align::VALTYPE type = get_type(M, swn);
                if (type == VALTYPE::NA && swn.s == 0)
                    GlobalTrack(Mval, M, swn, type, globalqueue, localqueue, valuequeue, id, dot_sources);
                else
                {
                    outputdot(Mval, M, swn, type);
                    switch (type)
                    {
                    case VALTYPE::Q:
                    {
                        SIZETYPE ib = 1;
                        for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
                        {
                            Node &node = node_map[nd];
                            if (!node.is_target)
                                continue;
                            if (ib > std::numeric_limits<BITSTYPE>::max())
                                break;
                            if (bits[M] & ib)
                                arrangesource(valuequeue, &node.Avals[compsz_map[nd] - 1][O.size()], localqueue, id);
                            ib *= 2;
                        }
                        break;
                    }
                    case VALTYPE::ABAR:
                        break;
                    case VALTYPE::SOE:
                    {
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Gvals[swn.s][swn.w - 1]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SOF:
                    {
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-Omax-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Gvals[swn.s - 1][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SOG:
                    {
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Evals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Fvals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & 4)
                        {
                            boost::graph_traits<graph_t>::vertex_descriptor td = boost::source(*eis[swn.n], graph);
                            arrangesource(valuequeue, &(node_map[td].Avals[compsz_map[td] - 1][swn.w]), localqueue, id);
                        }
                        if (bits[M] & 8)
                            arrangesource(valuequeue, &vals[M-Omax-2], localqueue, id);
                        break;
                    }
                    case VALTYPE::SID:
                    {
                        boost::graph_traits<graph_t>::vertex_descriptor td = boost::source(*eis[swn.n], graph);
                        arrangesource(valuequeue, &(node_map[td].Avals[compsz_map[td] - 1][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIE:
                    {
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Gvals[swn.s][swn.w - 1]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIF0:
                    {
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-Omax-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].G0vals[swn.s - 1][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIG0:
                    {
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Evals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].F0vals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(4))
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Gvals[swn.s - 1][swn.w - 1]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIG:
                    {
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].G0vals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].Dvals[swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(4))
                        {
                            boost::graph_traits<graph_t>::vertex_descriptor td = boost::source(*eis[swn.n], graph);
                            arrangesource(valuequeue, &(node_map[td].Avals[compsz_map[td] - 1][swn.w]), localqueue, id);
                        }
                        break;
                    }
                    case VALTYPE::LID0: case VALTYPE::SID0:
                        arrangesource(valuequeue, &(node_map[boost::source(*eis[swn.n], graph)].Avals[swn.s][swn.w]), localqueue, id);
                        break;
                    case VALTYPE::SIDX:
                    {
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &(node_map[boost::source(*eis[swn.n], graph)].Avals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(edge_map[*eis[swn.n]].D0vals[swn.s][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::NA:
                    {
                        boost::graph_traits<graph_t>::vertex_descriptor nd = Dot::DOTTYPE2NODEIDX(swn.n);
                        if (bits[M] & int8_t(1))
                            arrangesource(valuequeue, &node_map[nd].Avals[0][swn.w], localqueue, id);
                        SIZETYPE ib = 1;
                        boost::graph_traits<graph_t>::in_edge_iterator inei, inei_end;
                        for (boost::tie(inei, inei_end) = boost::in_edges(nd, graph); inei != inei_end; ++inei)
                        {
                            Edge &edge = edge_map[*inei];
                            if (comp_map[boost::source(*inei, graph)] != comp_map[nd])
                                continue;
                            ib *= 2;
                            if (ib > std::numeric_limits<BITSTYPE>::max())
                                break;
                            if (bits[M] & int8_t(ib))
                            {
                                SCORETYPE *source;
                                if (file2short.count(edge.name))
                                {
                                    if (Mval == edge.D0vals[swn.s - 1][swn.w] + edge.gfm) // this works because D0/DX[l-1] are only tracked by A[l]
                                        source = &(edge.D0vals[swn.s - 1][swn.w]);
                                    else
                                        source = &(edge.DXvals[swn.s - 1][swn.w]);
                                }
                                else
                                    source = &(edge.D0vals[swn.s - 1][swn.w]);
                                arrangesource(valuequeue, source, localqueue, id);
                            }
                        }
                        break;
                    }
                    case VALTYPE::NB:
                    {
                        boost::graph_traits<graph_t>::vertex_descriptor nd = Dot::DOTTYPE2NODEIDX(swn.n);
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &vals[M-1], localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(node_map[nd].Avals[compsz_map[nd] - 1][swn.w - 1]), localqueue, id);
                        break;
                    }
                    }
                }
            }
            else
            {
                Dot *M = globalqueue.front();
                globalqueue.pop();
                SCORETYPE Mval = valuequeue.front();
                valuequeue.pop();
                outputdot(Mval, M);
                for (SIZETYPE i = 0; i < M->s_sz; ++i)
                {
                    if (dot_sources[M->fs + i])
                        arrangesource(valuequeue, dot_sources[M->fs + i], globalqueue, id);
                    else
                    {
                        boost::graph_traits<graph_t>::vertex_descriptor td = boost::source(*eis[M->n], graph);
                        arrangesource(valuequeue, &node_map[td].Avals[compsz_map[td] - 1][M->w], localqueue, id);
                    }
                }
            }
        }
        if (id > max_id)
            max_id = id;

        for (boost::graph_traits<graph_t>::vertex_descriptor nd = 0; nd < boost::num_vertices(graph); ++nd)
            for (SIZETYPE w = 0; w <= O.size(); ++w)
            {
                node_map[nd].AdeltaDot[w].clear();
                node_map[nd].AdeltaGlobal[w].clear();
            }
    }

    void GlobalTrack(SCORETYPE Mval, SIZETYPE M, SWN &swn, Align::VALTYPE type, std::queue<Dot *> &globalqueue, std::queue<SIZETYPE> &localqueue, std::queue<SCORETYPE> &valuequeue, IDTYPE &id, std::deque<Dot *> &dot_sources)
    {
        boost::graph_traits<graph_t>::vertex_descriptor nd = Dot::DOTTYPE2NODEIDX(swn.n);
        Node &node = node_map[nd];
        outputdot(Mval, M, swn, type);
        for (SCORETYPE *source : node.AdeltaDot[swn.w])
            arrangesource(valuequeue, source, localqueue, id);

        for (Node::GlobalSuffix &globalsuffix : node.AdeltaGlobal[swn.w])
        {
            boost::graph_traits<graph_t>::edge_iterator ei = globalsuffix.ei;
            Edge &edge = edge_map[*ei];
            boost::graph_traits<graph_t>::vertex_descriptor hd = boost::target(*ei, graph), td = boost::source(*ei, graph);
            SCORETYPE *Avals = edge.tailAvals;
            TrackTree &tracktree = edge.tracktree;
            std::deque<TrackTree::TrackNode> &tracknodes = tracktree.tracknodes;
            RankVec &rankvec = file2rankvec[edge.name];
            std::ifstream &SAfin = file2SA[edge.name];
            SIZETYPE start;
            SAfin.seekg(globalsuffix.start * sizeof(SIZETYPE));
            SAfin.read((char*)&start, sizeof(SIZETYPE));
            std::ifstream &REVREFfin = file2long[edge.name];
            NUCTYPE revref[globalsuffix.lambda];
            REVREFfin.seekg(start);
            REVREFfin.read((char*)revref, globalsuffix.lambda);

            tracktree.pidx = -1;
            tracktree.idx = 0;
            tracktree.lambda = 0;
            if (comp_map[hd] != comp_map[td])
            {
                for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w == 0)
                    {
                        tracktree.dotE(w).val = -inf;
                        tracktree.dotE(w).s_sz = 0;
                    }
                    else
                        source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge.ue, tracktree.dotG(w - 1).val + edge.ve}, dot_sources);
                    tracktree.dotF(w).val = -inf;
                    tracktree.dotF(w).s_sz = 0;
                    source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w), NULL}, {tracktree.dotE(w).val, tracktree.dotF(w).val, Avals[w] + edge.T}, dot_sources);
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                for (SIZETYPE i = globalsuffix.lambda; i > 0;)
                {
                    --i;
                    NUCTYPE c = revref[i];
                    tracktree.goto_child(c, rankvec.bwt_sz - 1 - start - i);
                    for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                    {
                        if (w == 0)
                        {
                            tracktree.dotE(w).val = -inf;
                            tracktree.dotE(w).s_sz = 0;
                        }
                        else
                            source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge.ue, tracktree.dotG(w - 1).val + edge.ve}, dot_sources);
                        source_max(tracktree.dotF(w), {&tracktree.dotPF(w), &tracktree.dotPG(w)}, {tracktree.dotPF(w).val + edge.uf, tracktree.dotPG(w).val + edge.vf}, dot_sources);
                        if (w == 0)
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w)}, {tracktree.dotE(w).val, tracktree.dotF(w).val}, dot_sources);
                        else
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w), &tracktree.dotPG(w - 1)}, {tracktree.dotE(w).val, tracktree.dotF(w).val, tracktree.dotPG(w - 1).val + edge.gamma[c][O[w - 1]]}, dot_sources);
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                }
            }
            else
            {
                for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w > 0)
                        source_max(tracktree.dotG(w - 1), {&tracktree.dotE(w - 1), NULL}, {tracktree.dotE(w - 1).val, Avals[w - 1] + edge.T}, dot_sources);
                    if (w == 0)
                    {
                        tracktree.dotE(w).val = -inf;
                        tracktree.dotE(w).s_sz = 0;
                    }
                    else
                        source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge.ue, tracktree.dotG(w - 1).val + edge.ve}, dot_sources);
                    tracktree.dotF(w).val = -inf;
                    tracktree.dotF(w).s_sz = 0;
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                for (SIZETYPE i = globalsuffix.lambda; i > 0; )
                {
                    --i;
                    NUCTYPE c = revref[i];
                    tracktree.goto_child(c, rankvec.bwt_sz - 1 - start - i);
                    for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                    {
                        if (w == 0)
                        {
                            tracktree.dotE(w).val = -inf;
                            tracktree.dotE(w).s_sz = 0;
                        }
                        else
                            source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge.ue, tracktree.dotG(w - 1).val + edge.ve}, dot_sources);
                        if (tracktree.pidx == 0)
                            source_max(tracktree.dotF(w), {&tracktree.dotPF(w), &tracktree.dotPE(w)}, {tracktree.dotPF(w).val + edge.uf, tracktree.dotPE(w).val + edge.vf}, dot_sources);
                        else
                            source_max(tracktree.dotF(w), {&tracktree.dotPF(w), &tracktree.dotPG(w)}, {tracktree.dotPF(w).val + edge.uf, tracktree.dotPG(w).val + edge.vf}, dot_sources);
                        if (w == 0)
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w)}, {tracktree.dotE(w).val, tracktree.dotF(w).val}, dot_sources);
                        else
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w), &tracktree.dotPG(w - 1)}, {tracktree.dotE(w).val, tracktree.dotF(w).val, tracktree.dotPG(w - 1).val + edge.gamma[c][O[w - 1]]}, dot_sources);
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                }
            }
            arrangesource(valuequeue, &tracktree.dotG(swn.w), globalqueue, id);
        }
    }

    void source_max(Dot &dot, std::initializer_list<Dot *> adrs, std::initializer_list<SCORETYPE> vals, std::deque<Dot *> &dot_sources)
    {
        dot.val = -inf;
        for (std::initializer_list<SCORETYPE>::iterator itv = vals.begin(); itv != vals.end(); ++itv)
            if (*itv > dot.val)
                dot.val = *itv;
        dot.fs = dot_sources.size();
        std::initializer_list<Dot *>::iterator ita = adrs.begin();
        for (std::initializer_list<SCORETYPE>::iterator itv = vals.begin(); itv != vals.end(); ++itv, ++ita)
            if (*itv == dot.val)
                dot_sources.emplace_back(*ita);
        dot.s_sz = dot_sources.size() - dot.fs;
    }
};

#endif
