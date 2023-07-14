#ifndef ALIGN_H
#define ALIGN_H

#include <climits>
#include <new>
#include <initializer_list>
#include <iomanip>
#include <cmath>
#include <set>
#include <typeinfo>
#include "Graph.h"

struct TrackTree
{
    enum enumEFG
    {
        enumE,
        enumF,
        enumG
    };

    QUERYSIZE W;
    std::deque<Dot> dots;
    int64_t idx, pidx;
    SIZETYPE shiftE, shiftF, shiftG, shiftPE, shiftPF, shiftPG;

    void setidx(int64_t idx_)
    {
        if (idx_ == 0)
            pidx = -1;
        else
        {
            pidx = idx;
            shiftPE = shiftE;
            shiftPF = shiftF;
            shiftPG = shiftG;
        }
        idx = idx_;
        shiftE = 3 * idx * (W + 1);
        shiftF = shiftE + (W + 1);
        shiftG = shiftF + (W + 1);
    }
    struct TrackNode
    {
        int64_t cidxs[5];
        QUERYSIZE tau = 0;

        TrackNode(int64_t cidx_)
        {
            for (SIZETYPE i = 0; i < 5; ++i)
                cidxs[i] = cidx_;
        }
    };

    std::deque<TrackNode> tracknodes;

    void emplace_back(int64_t cidx_, DOTTYPE n_, SIZETYPE s_, QUERYSIZE lambda_)
    {
        tracknodes.emplace_back(cidx_);
        dots.resize(tracknodes.size() * (enumEFG::enumG + 1) * (W + 1));
        SIZETYPE k = (tracknodes.size() - 1) * (enumEFG::enumG + 1) * (W + 1);
        for (enumEFG t : {enumEFG::enumE, enumEFG::enumF, enumEFG::enumG})
            for (SIZETYPE w = 0; w <= W; ++w, ++k)
            {
                dots[k].n = n_;
                dots[k].s = s_;
                dots[k].w = w;
                dots[k].lambda = lambda_;
            }
    }

    void clear()
    {
        dots.clear();
        tracknodes.clear();
    }
};

template <typename T>
struct MonoDeque
{
    std::deque<T> data;
    SIZETYPE offset = 0;

    void emplace_back(T t)
    {
        if (offset < data.size())
            data[offset] = t;
        else
            data.emplace_back(t);
        ++offset;
    }

    T &operator[](SIZETYPE idx)
    {
        return data[idx];
    }

    void clear()
    {
        offset = 0;
    }
};

struct Align : Graph
{
    std::map<Edge *, TrackTree> long2tracktree;

    struct CrossGlobalData
    {
        struct Doo
        {
            QUERYSIZE w;
            SCORETYPE F, G;

            Doo(QUERYSIZE w_, SCORETYPE F_, SCORETYPE G_)
                : w(w_), F(F_), G(G_)
            {
            }
        };

        std::deque<Doo> FG;

        struct SimNode
        {
            NUCTYPE c;
            SIZETYPE sr1, sr2;
            SIZETYPE shiftFG;

            SimNode(NUCTYPE c_, SIZETYPE sr1_, SIZETYPE sr2_, SIZETYPE shiftFG_)
                : c(c_), sr1(sr1_), sr2(sr2_), shiftFG(shiftFG_)
            {
            }
        };

        std::deque<SimNode> simnodes;

        void emplace_back(NUCTYPE c_, SIZETYPE sr1_, SIZETYPE sr2_)
        {
            simnodes.emplace_back(c_, sr1_, sr2_, FG.size());
        }

        void pop_back()
        {
            if (FG.size() > simnodes.back().shiftFG)
                FG.erase(FG.begin() + simnodes.back().shiftFG, FG.end());
            simnodes.pop_back();
        }

        struct Do
        {
            QUERYSIZE w;
            SCORETYPE E;
        };

        std::unique_ptr<Do []> E0, E;
        std::unique_ptr<SCORETYPE []> Ethres;
    };

    std::mutex &mtx;
    std::ifstream &fin;
    std::vector<NUCTYPE> O;
    std::string Oname;
    std::ofstream fout;
    IDTYPE max_id = 0;
    QUERYSIZE Omax;
    std::vector<bool> Qbools;
    enum VALTYPE{Q, ABAR, LID0, SID0, SIDX, OTHERS};

    CrossGlobalData crossglobaldata;

    MonoDeque<Dot *> dot_sources;

    void apply_memory()
    {
        trn = 0;
        tnn = 1 + nodes.size();
        tsn = 0;

        std::vector<SIZETYPE> Asos;
        for (Node &node : nodes)
        {
            trn += node.scc_sz + 1;
            tnn += (Omax + 1) * (node.scc_sz + 1);
            Asos.push_back(1);
            for (EdgeLocalCircuit &local_circuit : local_circuits)
                if (local_circuit.head==&node)
                    Asos.back() += 2;
            for (EdgeGlobalCircuit &global_circuit : global_circuits)
                if (global_circuit.head==&node)
                    Asos.back() += 1;   
            tsn += node.get_tsn(Asos.back(), Omax);
        }
        for (EdgeLocalCross &edge : local_crosses)
        {
            SHORTSIZE ref_sz = file2short[edge.name].second;
            trn += edge.get_trn(ref_sz);
            tnn += edge.get_tnn(Omax, ref_sz);
            tsn += edge.get_tsn(Omax, ref_sz);
        }
        for (EdgeLocalCircuit &edge : local_circuits)
        {
            SHORTSIZE ref_sz = file2short[edge.name].second;
            trn += edge.get_trn(ref_sz);
            tnn += edge.get_tnn(Omax, ref_sz);
            tsn += edge.get_tsn(Omax, ref_sz);
        }
        sources.reset(new SCORETYPE *[tsn]);
        sourcess.reset(new SCORETYPE **[tnn-1-nodes.size()]);
        s_szs.reset(new SOURCESIZE[tnn-1-nodes.size()]);
        for (EdgeGlobalCircuit &edge : global_circuits)
        {
            trn += edge.get_trn();
            tnn += edge.get_tnn(Omax);
        }
        vals.reset(new SCORETYPE[tnn]);
        ids.reset(new IDTYPE[tnn]);
        ss.reset(new SHORTSIZE[trn]);
        ns.reset(new DOTTYPE[trn]);

        SCORETYPE *fpval = vals.get(), *rpval = fpval + tnn;
        pQval = --rpval;
        SCORETYPE **fpsource = sources.get();
        SCORETYPE ***fpsources = sourcess.get();
        SOURCESIZE *fps_sz = s_szs.get();
        IDTYPE *fpid = ids.get(), *rpid = fpid + tnn;
        pQid = --rpid;
        SHORTSIZE *fps = ss.get();
        DOTTYPE *fpn = ns.get();
        for (SIZETYPE i=0; i<nodes.size(); ++i)
            nodes[i].apply_memory(Asos[i], Omax, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn, rpval, rpid);
        for (EdgeLocalCross &edge : local_crosses)
            edge.apply_memory(Omax, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn, file2short[edge.name].second);
        for (EdgeLocalCircuit &edge : local_circuits)
            edge.apply_memory(Omax, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn, file2short[edge.name].second);
        for (EdgeGlobalCircuit &edge : global_circuits)
            edge.apply_memory(Omax, fpval, fpid, fps, fpn);
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

    Align(std::mutex &mtx_, std::ifstream &fin_, boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short, std::map<std::string, RankVec> &file2rankvec, std::string mgfile, QUERYSIZE Omax_)
        : mtx(mtx_), fin(fin_), Graph(vm, file2short, file2rankvec), fout(mgfile, std::ifstream::binary), Omax(Omax_)
    {
        for (Node &node : nodes)
            if (node.is_target)
                Qbools.push_back(false);
        

        apply_memory();

        crossglobaldata.E0.reset(new CrossGlobalData::Do[Omax + 1]);
        crossglobaldata.E.reset(new CrossGlobalData::Do[Omax + 1]);
        crossglobaldata.Ethres.reset(new SCORETYPE[Omax + 1]);
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
        for (SIZETYPE i = 0; i < nodes.size(); ++i)
        {
            for (SIZETYPE w = 0; w <= O.size(); ++w)
                nodes[i].Avals[0][w] = -inf;
            if (nodes[i].is_root)
                *nodes[i].pAbarval = 0;
            else
                *nodes[i].pAbarval = -inf;
        }
        for (SIZETYPE scc_id = 0, scc_sz = 0, cum_scc_sz = 0; cum_scc_sz < nodes.size(); ++scc_id, cum_scc_sz += scc_sz, scc_sz = 0)
        {
            for (Node &node : nodes)
                if (node.scc_id == scc_id)
                    ++scc_sz;

            for (SIZETYPE w = 0; w <= O.size(); ++w)
            {
                for (Node &node : nodes)
                {
                    if (node.scc_id != scc_id)
                        continue;

                    if (w == 0)
                        source_max(&node.Bvals[w], {}, {});
                    else
                        source_max(&node.Bvals[w], {&node.Bvals[w - 1], &node.Avals[node.scc_sz - 1][w - 1]}, {node.Bvals[w - 1] + node.ue, node.Avals[node.scc_sz - 1][w - 1] + node.ve});

                    node.updateA0(w, node.Bvals[w], &node.Bvals[w]);
                    if (w == 0)
                        node.updateA0(w, node.pAbarval[0], node.pAbarval);
                }
                for (Edge *edge : edges)
                    if (edge->head->scc_id == scc_id && edge->tail->scc_id == scc_id && file2short.count(edge->name))
                        CircuitIteration(w, (EdgeLocalCircuit*)edge);
                CircuitIterationGlobal(w, scc_id);
                for (SIZETYPE l = 1; l < scc_sz; ++l)
                {
                    for (Edge *edge : edges)
                    {
                        if (edge->head->scc_id != scc_id || edge->tail->scc_id != scc_id || !file2long.count(edge->name))
                            continue;
                        EdgeGlobalCircuit *circuit = (EdgeGlobalCircuit *)edge;
                        circuit->D0vals[l - 1][w] = circuit->tail->Avals[l - 1][w] + circuit->T;
                    }
                    for (Edge *edge : edges)
                    {
                        if (edge->head->scc_id != scc_id || edge->tail->scc_id != scc_id || !file2short.count(edge->name))
                            continue;
                        EdgeLocalCircuit *circuit = (EdgeLocalCircuit *)edge;
                        SHORTSIZE ref_sz = file2short[circuit->name].second;
                        circuit->D0vals[l - 1][w] = circuit->tail->Avals[l - 1][w] + circuit->T;
                        circuit->DXvals[l - 1][w] = std::max(circuit->tail->Avals[l - 1][w] + circuit->gfpT[ref_sz], circuit->D0vals[l - 1][w] + circuit->gf[ref_sz]);
                        if (circuit->DXvals[l - 1][w] == circuit->tail->Avals[l - 1][w] + circuit->gfpT[ref_sz])
                            circuit->DXbits[(l-1)*(Omax+1)+w] |= uint8_t(1);
                        if (circuit->DXvals[l - 1][w] == circuit->D0vals[l - 1][w] + circuit->gf[ref_sz])
                            circuit->DXbits[(l-1)*(Omax+1)+w] |= uint8_t(2);
                    }
                    for (Node &node : nodes)
                    {
                        if (node.scc_id != scc_id)
                            continue;
                        
                        node.Avals[l][w] = node.Avals[0][w];
                        for (Edge *edge : edges)
                        {
                            if (edge->head != &node || edge->tail->scc_id != scc_id || !file2short.count(edge->name))
                                continue;
                            EdgeLocalCircuit *circuit = (EdgeLocalCircuit *)edge;
                            node.Avals[l][w] = std::max({node.Avals[l][w], circuit->D0vals[l - 1][w] + circuit->gfm, circuit->DXvals[l - 1][w]});
                        }
                        for (Edge *edge : edges)
                        {
                            if (edge->head != &node || edge->tail->scc_id != scc_id || !file2long.count(edge->name))
                                continue;
                            EdgeGlobalCircuit *circuit = (EdgeGlobalCircuit *)edge;
                            node.Avals[l][w] = std::max(node.Avals[l][w], circuit->D0vals[l - 1][w]);
                        }
                        SIZETYPE idx = &node.Avals[l][w] - vals.get();
                        s_szs[idx] = 0;
                        if (node.Avals[0][w] == node.Avals[l][w])
                            sourcess[idx][s_szs[idx]++] = &node.Avals[0][w];
                        for (Edge *edge : edges)
                        {
                            if (edge->head != &node || edge->tail->scc_id != scc_id || !file2short.count(edge->name))
                                continue;
                            EdgeLocalCircuit *circuit = (EdgeLocalCircuit *)edge;
                            if (circuit->D0vals[l - 1][w] + circuit->gfm == node.Avals[l][w])
                                sourcess[idx][s_szs[idx]++] = &circuit->D0vals[l - 1][w];
                            if (circuit->DXvals[l - 1][w] == node.Avals[l][w])
                                sourcess[idx][s_szs[idx]++] = &circuit->DXvals[l - 1][w];
                        }
                        for (Edge *edge : edges)
                        {
                            if (edge->head != &node || edge->tail->scc_id !=scc_id || !file2long.count(edge->name))
                                continue;
                            EdgeGlobalCircuit *circuit = (EdgeGlobalCircuit *)edge;
                            if (circuit->D0vals[l - 1][w] == node.Avals[l][w])
                                sourcess[idx][s_szs[idx]++] = &circuit->D0vals[l - 1][w];
                        }
                    }
                }
            }
            for (Edge *edge : edges)
                if (edge->tail->scc_id == scc_id && edge->head->scc_id != scc_id && file2short.count(edge->name))
                    CrossIteration((EdgeLocalCross *)edge);
            for (Edge *edge : edges)
                if (edge->tail->scc_id == scc_id && edge->head->scc_id != scc_id && file2long.count(edge->name))
                    CrossIterationGlobal(edge);
        }
        *pQval = -inf;
        for (Node &node : nodes)
            if (node.is_target)
                *pQval = std::max(*pQval, node.Avals[node.scc_sz - 1][O.size()]);
        for (SIZETYPE i = 0, j = 0; i < nodes.size(); ++i)
            if (nodes[i].is_target && nodes[i].Avals[nodes[i].scc_sz - 1][O.size()] == *pQval)
                Qbools[j++] = true;
    }

    void CrossIteration(EdgeLocalCross *edge)
    {
        NUCTYPE *ref = file2short[edge->name].first.get();
        SHORTSIZE ref_sz = file2short[edge->name].second;
        SCORETYPE *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        SCORETYPE **Evals = edge->Evals.get(), **Fvals = edge->Fvals.get(), **Gvals = edge->Gvals.get();
        for (SIZETYPE s = 0; s <= ref_sz; ++s)
            for (SIZETYPE w = 0; w <= O.size(); ++w)
            {
                if (w == 0)
                    source_max(&Evals[s][w], {}, {});
                else
                    source_max(&Evals[s][w], {&Evals[s][w - 1], &Gvals[s][w - 1]}, {Evals[s][w - 1] + edge->ue, Gvals[s][w - 1] + edge->ve});
                if (s == 0)
                    source_max(&Fvals[s][w], {}, {});
                else
                    source_max(&Fvals[s][w], {&Fvals[s - 1][w], &Gvals[s - 1][w]}, {Fvals[s - 1][w] + edge->uf, Gvals[s - 1][w] + edge->vf});
                if (s == 0 || w == 0)
                    source_max(&Gvals[s][w], {&Fvals[s][w], &Evals[s][w], &Avals[w]}, {Fvals[s][w], Evals[s][w], Avals[w] + edge->gfpT[s]});
                else
                    source_max(&Gvals[s][w], {&Fvals[s][w], &Evals[s][w], &Gvals[s - 1][w - 1], &Avals[w]}, {Fvals[s][w], Evals[s][w], Gvals[s - 1][w - 1] + edge->gamma[ref[s - 1]][O[w - 1]], Avals[w] + edge->gfpT[s]});
                edge->head->updateA0(w, Gvals[s][w], &Gvals[s][w]);
            }
    }

    void CrossIterationGlobal(Edge *edge)
    {
        SCORETYPE *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        std::deque<CrossGlobalData::SimNode> &simnodes = crossglobaldata.simnodes;
        std::unique_ptr<CrossGlobalData::Do[]> &E0 = crossglobaldata.E0, &E = crossglobaldata.E;
        std::unique_ptr<SCORETYPE []> &Ethres = crossglobaldata.Ethres;
        std::deque<CrossGlobalData::Doo> &FG = crossglobaldata.FG;
        RankVec &rankvec = file2rankvec[edge->name]; 

        SIZETYPE sr1 = 0, sr2 = rankvec.bwt_sz;
        crossglobaldata.emplace_back(0, sr1, sr2); // the root SimNode has no base, we simply set it to 0, which does not mean that it has base #(0) 
        for (SIZETYPE w = 0; w <= O.size(); ++w)
        {
            SCORETYPE tgw = (O.size() > w ? (O.size() - w - 1) * edge->tail->ue + edge->tail->ve : 0);
            E0[w].w = w;
            if (w > 0)
                E0[w].E = std::max(E0[w - 1].E + edge->ue, FG[w - 1].G + edge->ve);
            else
                E0[w].E = -inf;
            SCORETYPE tmp = std::max(E0[w].E, Avals[w] + edge->T);
            if (E0[w].E + edge->ue < tmp + edge->tail->ve && E0[w].E + (O.size() - w) * edge->ue < tmp + tgw)
                E0[w].E = -inf;
            FG.emplace_back(w, -inf, std::max(E0[w].E, Avals[w] + edge->T));
            Ethres[w] = std::max(E0[w].E, FG[w].G + std::min({SCORETYPE(0), edge->tail->ve - edge->ue, SCORETYPE(tgw - (O.size() - w) * edge->ue)}));
            edge->head->updateA0(w, FG[w].G, edge, simnodes[0].sr1, 0);
        }

        while (true)
        {
            do
            {
                NUCTYPE c;
                if (simnodes.back().shiftFG < FG.size())
                    c = 1;
                else
                {
                    do
                    {
                        c = simnodes.back().c;
                        crossglobaldata.pop_back();
                    } while (c == 6 && !simnodes.empty());
                    if (simnodes.empty())
                        break;
                }
                do
                {
                    ++c;
                    sr1 = rankvec.C[c] + rankvec.rank(c, simnodes.back().sr1);
                    sr2 = rankvec.C[c] + rankvec.rank(c, simnodes.back().sr2);
                } while (sr1 >= sr2 && c < 6);
                if (sr1 < sr2)
                {
                    crossglobaldata.emplace_back(c, sr1, sr2);
                    break;
                }
                else if (FG.size() > simnodes.back().shiftFG)
                    FG.erase(FG.begin() + simnodes.back().shiftFG, FG.end());
            } while (true);
            if (simnodes.empty())
                break;

            SIZETYPE idx = simnodes[simnodes.size() - 2].shiftFG;
            for (SIZETYPE e = 0, w = FG[idx].w;; ++e)
            {
                E[e].w = w;
                E[e].E = std::max((e > 0 && E[e - 1].w == w - 1) ? E[e - 1].E + edge->ue : -inf, (simnodes.back().shiftFG < FG.size() && FG.back().w == w - 1) ? FG.back().G + edge->ve : -inf);
                if (E[e].E < Ethres[w])
                    E[e].E = -inf;

                SCORETYPE F_val = idx < simnodes.back().shiftFG && FG[idx].w == w ? std::max(FG[idx].F + edge->uf, FG[idx].G + edge->vf) : -inf;
                SCORETYPE G_val = idx > simnodes[simnodes.size() - 2].shiftFG && FG[idx - 1].w == w - 1 ? FG[idx - 1].G + edge->gamma[simnodes.back().c][O[w - 1]] : -inf;
                G_val = std::max({G_val, E[e].E, F_val});
                if (G_val >= FG[w].G)
                {
                    FG.emplace_back(w, F_val, G_val);
                    edge->head->updateA0(w, G_val, edge, simnodes.back().sr1, simnodes.size() - 1);
                }
                if (idx < simnodes.back().shiftFG && FG[idx].w <= w)
                    ++idx;
                if (w < O.size() && (FG[idx - 1].w == w || (simnodes.back().shiftFG < FG.size() && FG.back().w == w) || E[e].E > -inf))
                    ++w;
                else if (idx < simnodes.back().shiftFG)
                    w = FG[idx].w;
                else
                    break;
            }
        }
    }

    void CircuitIteration(QUERYSIZE w, EdgeLocalCircuit *edge)
    {
        NUCTYPE *ref = file2short[edge->name].first.get();
        SHORTSIZE ref_sz = file2short[edge->name].second;
        SCORETYPE *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        SCORETYPE *Dvals = edge->Dvals;
        SCORETYPE **Evals = edge->Evals.get(), **F0vals = edge->F0vals.get(), **G0vals = edge->G0vals.get(), **Gvals = edge->Gvals.get();

        if (w > 0)
            source_max(&Dvals[w - 1], {&Avals[w - 1]}, {Avals[w - 1] + edge->T});
        for (SIZETYPE s = 0; s <= ref_sz; ++s)
        {
            if (w > 0)
            {
                if (s > 0)
                    source_max(&Gvals[s][w - 1], {&G0vals[s][w - 1], &Dvals[w - 1], &Avals[w - 1]}, {G0vals[s][w - 1], Dvals[w - 1] + edge->gf[s], Avals[w - 1] + edge->gfpT[s]});
                else
                    source_max(&Gvals[s][w - 1], {&G0vals[s][w - 1], &Dvals[w - 1]}, {G0vals[s][w - 1], Dvals[w - 1]});
            }
            if (w == 0)
                source_max(&Evals[s][w], {}, {});
            else
                source_max(&Evals[s][w], {&Evals[s][w - 1], &Gvals[s][w - 1]}, {Evals[s][w - 1] + edge->ue, Gvals[s][w - 1] + edge->ve});
            if (s == 0)
                source_max(&F0vals[s][w], {}, {});
            else
                source_max(&F0vals[s][w], {&F0vals[s - 1][w], &G0vals[s - 1][w]}, {F0vals[s - 1][w] + edge->uf, G0vals[s - 1][w] + edge->vf});
            if (s > 0 && w > 0)
                source_max(&G0vals[s][w], {&Evals[s][w], &F0vals[s][w], &Gvals[s - 1][w - 1]}, {Evals[s][w], F0vals[s][w], Gvals[s - 1][w - 1] + edge->gamma[ref[s - 1]][O[w - 1]]});
            else
                source_max(&G0vals[s][w], {&Evals[s][w], &F0vals[s][w]}, {Evals[s][w], F0vals[s][w]});
            edge->head->updateA0(w, G0vals[s][w], &G0vals[s][w]);
        }
    }

    void CircuitIterationGlobal(QUERYSIZE w, SIZETYPE scc_id)
    {
        if (w == 0)
            return;

        std::vector<EdgeGlobalCircuit *> circuits;
        for (Edge *edge : edges)
            if (edge->head->scc_id == scc_id && edge->tail->scc_id == scc_id && file2long.count(edge->name))
                circuits.push_back((EdgeGlobalCircuit *)edge);
        std::vector<SCORETYPE> tgws(circuits.size());
        std::vector<SNC *> jumps(circuits.size());
        std::vector<RankVec *> prankvecs;
        for (SIZETYPE i = 0; i < circuits.size(); ++i)
            prankvecs.push_back(&file2rankvec[circuits[i]->name]);

        for (SIZETYPE i = 0; i < circuits.size(); ++i)
        {
            if (w == 1)
            {
                circuits[i]->sncs.emplace_back(-inf, -inf, -inf, -inf, 0, 5, 0, 0, prankvecs[i]->bwt_sz, 1, (SNC *)NULL, (SNC *)NULL); // the root SNC has no base, we simply set it to 0, which does not mean that it has base #(0)
                for (NUCTYPE c = 2; c <= 6; ++c)
                {
                    SIZETYPE sr1 = prankvecs[i]->C[c] + prankvecs[i]->rank(c, circuits[i]->sncs.front().sr1);
                    SIZETYPE sr2 = prankvecs[i]->C[c] + prankvecs[i]->rank(c, circuits[i]->sncs.front().sr2);
                    if (sr1 < sr2)
                    {
                        circuits[i]->sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, 1, sr1, sr2, 0, &(circuits[i]->sncs.front()), circuits[i]->sncs.front().jump);
                        circuits[i]->sncs.front().jump = &(circuits[i]->sncs.back());
                    }
                }
            }
            else
            {
                for (SNC *prejump = &(circuits[i]->sncs.front()), *jump = prejump->jump; jump; jump = prejump->jump)
                {
                    jump->hatG = jump->G0;
                    if (jump->itp != &(circuits[i]->sncs.front()) && jump->itp->hatG == -inf && jump->hatG == -inf && jump->E == -inf)
                        prejump->jump = jump->jump;
                    else
                        prejump = jump;
                }
            }

            circuits[i]->sncs.front().hatG = std::max(circuits[i]->sncs.front().G0, circuits[i]->tail->Avals[circuits[i]->tail->scc_sz - 1][w - 1] + circuits[i]->T);
            circuits[i]->sncs.front().E = std::max(circuits[i]->sncs.front().hatG + circuits[i]->ve, circuits[i]->sncs.front().E + circuits[i]->ue);
            tgws[i] = (O.size() > w ? (O.size() - w - 1) * circuits[i]->tail->ue + circuits[i]->tail->ve : 0);
            circuits[i]->sncs.front().F0 = -inf;
            circuits[i]->sncs.front().G0 = circuits[i]->sncs.front().E;
            circuits[i]->head->updateA0(w, circuits[i]->sncs.front().G0, circuits[i], circuits[i]->sncs.front().sr1, circuits[i]->sncs.front().lambda);

            jumps[i] = circuits[i]->sncs.front().jump;
        }

        while (circuits.size() > 0)
            for (SIZETYPE i = 0; i < circuits.size();)
            {
                SCORETYPE Gthres = std::max(circuits[i]->sncs.front().E, circuits[i]->tail->Avals[0][w] + circuits[i]->T);
                SCORETYPE Ethres = std::max(circuits[i]->sncs.front().E, Gthres + std::min(SCORETYPE(0), std::min(circuits[i]->tail->ve - circuits[i]->ue, SCORETYPE(tgws[i] - (O.size() - w) * circuits[i]->ue))));

                if (jumps[i])
                {
                    jumps[i]->E = std::max(jumps[i]->hatG + circuits[i]->ve, jumps[i]->E + circuits[i]->ue);
                    if (jumps[i]->E < Ethres)
                        jumps[i]->E = -inf;
                    jumps[i]->F0 = std::max(jumps[i]->itp->G0 + circuits[i]->vf, jumps[i]->itp->F0 + circuits[i]->uf);
                    jumps[i]->G0 = std::max(std::max(jumps[i]->E, jumps[i]->F0), jumps[i]->itp->hatG + circuits[i]->gamma[jumps[i]->c][O[w - 1]]);
                    if (jumps[i]->G0 < Gthres)
                    {
                        jumps[i]->F0 = -inf;
                        jumps[i]->G0 = -inf;
                    }
                    else
                        circuits[i]->head->updateA0(w, jumps[i]->G0, circuits[i], jumps[i]->sr1, jumps[i]->lambda);
                    if (jumps[i]->G0 != -inf && jumps[i]->hatG == -inf)
                    {
                        if (jumps[i]->cid == 0)
                        {
                            jumps[i]->cid = circuits[i]->sncs.size();
                            for (NUCTYPE c = 2; c < RankVec::sigma; ++c)
                            {
                                SIZETYPE sr1 = prankvecs[i]->C[c] + prankvecs[i]->rank(c, jumps[i]->sr1);
                                SIZETYPE sr2 = prankvecs[i]->C[c] + prankvecs[i]->rank(c, jumps[i]->sr2);
                                if (sr1 < sr2)
                                {
                                    circuits[i]->sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, jumps[i]->lambda + 1, sr1, sr2, 0, jumps[i], (SNC *)NULL);
                                    ++jumps[i]->idl;
                                }
                            }
                        }
                        for (NUCTYPE idi = 0; idi < jumps[i]->idl; ++idi)
                            if (circuits[i]->sncs[jumps[i]->cid + idi].E == -inf && circuits[i]->sncs[jumps[i]->cid + idi].G0 == -inf)
                            {
                                circuits[i]->sncs[jumps[i]->cid + idi].jump = jumps[i]->jump;
                                jumps[i]->jump = &(circuits[i]->sncs[jumps[i]->cid + idi]);
                            }
                    }
                    jumps[i] = jumps[i]->jump;
                    
                    ++ i;
                }
                else
                {
                    if (w == O.size())
                        circuits[i]->sncs.clear();
                    circuits.erase(circuits.begin() + i);
                    tgws.erase(tgws.begin() + i);
                    jumps.erase(jumps.begin() + i);
                }
            }
    }

    Align::VALTYPE get_type(SIZETYPE M, SWN &swn)
    {
        if (swn.n==Dot::DotQ)
            return VALTYPE::Q;
        if (swn.n==Dot::DotAbar)
            return VALTYPE::ABAR;
        if (swn.n>=0)
        {
            if (file2long.count(edges[swn.n]->name))
                return VALTYPE::LID0;
            if (edges[swn.n]->tail->scc_id==edges[swn.n]->head->scc_id)
            {
                EdgeLocalCircuit *egcr = (EdgeLocalCircuit *) edges[swn.n];
                if (M >= egcr->D0vals[0]-vals.get())
                {
                    if (M < egcr->DXvals[0]-vals.get())
                        return VALTYPE::SID0;
                    return VALTYPE::SIDX;
                }
            }
        }
        return VALTYPE::OTHERS;
    }

    void outputdot(SCORETYPE Mval, SIZETYPE M, SWN &swn, Align::VALTYPE type)
    {
        fout.write((char *)&swn.n, sizeof(Dot::n));
        fout.write((char *)&swn.s, sizeof(Dot::s));
        fout.write((char *)&swn.w, sizeof(Dot::w));
        fout.write((char *)&Mval, sizeof(Dot::val));
        SOURCESIZE s_sz;
        if (type==VALTYPE::Q)
        {
            s_sz = 0;
            for (SIZETYPE i = 0; i < Qbools.size(); ++i)
                if (Qbools[i])
                    ++s_sz;
        }
        else
        {
            if (type==VALTYPE::ABAR)
                s_sz = 0;
            else
            {
                if (type==VALTYPE::LID0 || type==VALTYPE::SID0)
                    s_sz = 1;
                else
                {
                    if (type==VALTYPE::SIDX)
                    {
                        EdgeLocalCircuit *egcr = (EdgeLocalCircuit *) edges[swn.n];
                        s_sz = (egcr->DXbits[swn.s*(Omax+1)+swn.w] + 1) / 2;
                    }
                    else
                        s_sz = s_szs[M];
                }
            }
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
        std::map<Edge *, std::unique_ptr<SCORETYPE []>> edge2tailAvals;
        for (Edge *edge : edges)
        {
            edge2tailAvals[edge].reset(new SCORETYPE[O.size()+1]);
            for (SIZETYPE i=0; i<=O.size(); ++i)
                edge2tailAvals[edge].get()[i] = edge->tail->Avals[edge->tail->scc_sz - 1][i];
        }

        SIZETYPE Onsz = Oname.size();
        fout.write((char *)&Onsz, sizeof(Onsz));
        fout.write(Oname.data(), Onsz);
        std::queue<Dot *> globalqueue;
        std::queue<SIZETYPE> localqueue;
        std::queue<SCORETYPE> valuequeue;
        IDTYPE id = 0;
        *pQid = id;
        localqueue.push(pQval - vals.get());
        valuequeue.push(*pQval);
        *pQval = -inf;
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
                if (swn.n < Dot::DotQ && swn.s == 0)
                    GlobalTrack(edge2tailAvals, Mval, M, swn, type, globalqueue, localqueue, valuequeue, id);
                else
                {
                    outputdot(Mval, M, swn, type);
                    if (type==VALTYPE::Q)
                    {
                        for (SIZETYPE i = 0, j = 0; i < nodes.size(); ++i)
                            if (nodes[i].is_target && Qbools[j])
                            {
                                arrangesource(valuequeue, &(nodes[i].Avals[nodes[i].scc_sz - 1][O.size()]), localqueue, id);
                                Qbools[j++] = false;
                            }
                    }
                    else
                    {
                        if (type==VALTYPE::LID0 || type==VALTYPE::SID0)
                            arrangesource(valuequeue, &(edges[swn.n]->tail->Avals[swn.s][swn.w]), localqueue, id);
                        else
                        {
                            if (type==VALTYPE::SIDX)
                            {
                                EdgeLocalCircuit *egcr = (EdgeLocalCircuit *) edges[swn.n];
                                uint8_t &DXbit = egcr->DXbits[swn.s*(Omax+1)+swn.w];
                                if (DXbit & uint8_t(1))
                                    arrangesource(valuequeue, &(egcr->tail->Avals[swn.s][swn.w]), localqueue, id);
                                if (DXbit & uint8_t(2))
                                    arrangesource(valuequeue, &(egcr->D0vals[swn.s][swn.w]), localqueue, id);
                                DXbit = 0;
                            }
                            else if (type!=VALTYPE::ABAR)
                                for (SIZETYPE i = 0; i < s_szs[M]; ++i)
                                    arrangesource(valuequeue, sourcess[M][i], localqueue, id);
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
                        Node *tail = edges[M->n]->tail;
                        arrangesource(valuequeue, &tail->Avals[tail->scc_sz - 1][M->w], localqueue, id);
                    }
                }
            }
        }
        if (id > max_id)
            max_id = id;

        for (std::map<Edge *, TrackTree>::iterator it = long2tracktree.begin(); it != long2tracktree.end(); ++it)
            it->second.clear();
        for (Node &node : nodes)
            node.clearAdelta(O.size());
        dot_sources.clear();
    }

    void GlobalTrack(std::map<Edge *, std::unique_ptr<SCORETYPE []>> &edge2tailAvals, SCORETYPE Mval, SIZETYPE M, SWN &swn, Align::VALTYPE type, std::queue<Dot *> &globalqueue, std::queue<SIZETYPE> &localqueue, std::queue<SCORETYPE> &valuequeue, IDTYPE &id)
    {
        Node &node = nodes[Dot::nidx_trans(swn.n)];
        s_szs[M] = node.AdeltaGlobal[swn.w].size() + node.AdeltaDot[swn.w].size();
        outputdot(Mval, M, swn, type);
        for (SCORETYPE *source : node.AdeltaDot[swn.w])
            arrangesource(valuequeue, source, localqueue, id);

        for (Node::GlobalSuffix &globalsuffix : node.AdeltaGlobal[swn.w])
        {
            Edge *edge = globalsuffix.edge;
            SCORETYPE *Avals = edge2tailAvals[edge].get();
            TrackTree &tracktree = long2tracktree[edge];
            std::deque<TrackTree::TrackNode> &tracknodes = tracktree.tracknodes;
            std::deque<Dot> &dots = tracktree.dots;
            RankVec &rankvec = file2rankvec[edge->name];
            std::ifstream &SAfin = file2SA[edge->name];
            SIZETYPE start;
            SAfin.seekg(globalsuffix.start * sizeof(SIZETYPE));
            SAfin.read((char*)&start, sizeof(SIZETYPE));
            std::ifstream &REVREFfin = file2long[edge->name];
            std::unique_ptr<NUCTYPE []> revref(new NUCTYPE[globalsuffix.lambda]);
            REVREFfin.seekg(start);
            REVREFfin.read((char*)revref.get(), globalsuffix.lambda);

            if (tracktree.tracknodes.empty())
            {
                tracktree.W = O.size();
                tracktree.emplace_back(-1, edge->n, rankvec.bwt_sz - 1 - start - globalsuffix.lambda, 0); // minus
            }
            tracktree.setidx(0);
            if (edge->head->scc_id != edge->tail->scc_id)
            {
                for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w == 0)
                        source_max(dots[tracktree.shiftE + w], {}, {});
                    else
                        source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                    source_max(dots[tracktree.shiftF + w], {}, {});
                    source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], NULL}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, Avals[w] + edge->T});
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                for (SIZETYPE i = globalsuffix.lambda; i > 0;)
                {
                    --i;
                    NUCTYPE c = revref[i];
                    if (tracknodes[tracktree.idx].cidxs[c - 2] < 0)
                    {
                        tracktree.emplace_back(-1, edge->n, rankvec.bwt_sz - 1 - start - i, dots[tracktree.shiftE].lambda + 1);
                        tracknodes[tracktree.idx].cidxs[c - 2] = tracknodes.size() - 1;
                    }
                    tracktree.setidx(tracknodes[tracktree.idx].cidxs[c - 2]);
                    for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                    {
                        if (w == 0)
                            source_max(dots[tracktree.shiftE + w], {}, {});
                        else
                            source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                        source_max(dots[tracktree.shiftF + w], {&dots[tracktree.shiftPF + w], &dots[tracktree.shiftPG + w]}, {dots[tracktree.shiftPF + w].val + edge->uf, dots[tracktree.shiftPG + w].val + edge->vf});
                        if (w == 0)
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val});
                        else
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], &dots[tracktree.shiftPG + w - 1]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, dots[tracktree.shiftPG + w - 1].val + edge->gamma[c][O[w - 1]]});
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                }
            }
            else
            {
                for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w > 0)
                        source_max(dots[tracktree.shiftG + w - 1], {&dots[tracktree.shiftE + w - 1], NULL}, {dots[tracktree.shiftE + w - 1].val, Avals[w - 1] + edge->T});
                    if (w == 0)
                        source_max(dots[tracktree.shiftE + w], {}, {});
                    else
                        source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                    source_max(dots[tracktree.shiftF + w], {}, {});
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                for (SIZETYPE i = globalsuffix.lambda; i > 0; )
                {
                    --i;
                    NUCTYPE c = revref[i];
                    if (tracknodes[tracktree.idx].cidxs[c - 2] < 0)
                    {
                        tracktree.emplace_back(-1, edge->n, rankvec.bwt_sz - 1 - start - i, dots[tracktree.shiftE].lambda + 1); // minus
                        tracknodes[tracktree.idx].cidxs[c - 2] = tracknodes.size() - 1;
                    }
                    tracktree.setidx(tracknodes[tracktree.idx].cidxs[c - 2]);
                    for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                    {
                        if (w == 0)
                            source_max(dots[tracktree.shiftE + w], {}, {});
                        else
                            source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                        if (tracktree.pidx == 0)
                            source_max(dots[tracktree.shiftF + w], {&dots[tracktree.shiftPF + w], &dots[tracktree.shiftPE + w]}, {dots[tracktree.shiftPF + w].val + edge->uf, dots[tracktree.shiftPE + w].val + edge->vf});
                        else
                            source_max(dots[tracktree.shiftF + w], {&dots[tracktree.shiftPF + w], &dots[tracktree.shiftPG + w]}, {dots[tracktree.shiftPF + w].val + edge->uf, dots[tracktree.shiftPG + w].val + edge->vf});
                        if (w == 0)
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val});
                        else
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], &dots[tracktree.shiftPG + w - 1]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, dots[tracktree.shiftPG + w - 1].val + edge->gamma[c][O[w - 1]]});
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                }
            }
            arrangesource(valuequeue, &dots[tracktree.shiftG + swn.w], globalqueue, id);
        }
    }

    void source_max(SCORETYPE *val, std::initializer_list<SCORETYPE *> sources, std::initializer_list<SCORETYPE> source_vals)
    {
        SIZETYPE idx = val - vals.get();
        s_szs[idx] = 0;
        if (source_vals.begin() == source_vals.end())
            *val = -inf;
        else
        {
            *val = std::max(source_vals);
            std::initializer_list<SCORETYPE *>::iterator sources_it = sources.begin();
            for (auto source_vals_it = source_vals.begin(); source_vals_it != source_vals.end(); ++source_vals_it, ++sources_it)
                if (*source_vals_it == *val)
                    sourcess[idx][s_szs[idx]++] = *sources_it;
        }
    }

    void source_max(Dot &dot, std::initializer_list<Dot *> adrs, std::initializer_list<SCORETYPE> vals)
    {
        source_max(dot, adrs.begin(), vals.begin(), vals.end());
    }

    template <typename Ita, typename Itv>
    void source_max(Dot &dot, Ita ita_b, Itv itv_b, Itv itv_e)
    {
        dot.val = -inf;
        for (Itv itv = itv_b; itv != itv_e; ++itv)
            if (*itv > dot.val)
                dot.val = *itv;
        dot.fs = dot_sources.offset;
        for (Itv itv = itv_b; itv != itv_e; ++itv, ++ita_b)
            if (*itv == dot.val)
                dot_sources.emplace_back(*ita_b);
        dot.s_sz = dot_sources.offset - dot.fs;
    }
};

#endif
