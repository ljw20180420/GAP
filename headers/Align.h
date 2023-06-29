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

template <typename T>
struct MonoDeque
{
    std::deque<T> data;
    uint64_t offset = 0;

    void emplace_back(T t)
    {
        if (offset < data.size())
            data[offset] = t;
        else
            data.emplace_back(t);
        ++offset;
    }

    T &operator[](int64_t idx)
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
    struct CrossGlobalData
    {
        struct Doo
        {
            int w;
            double F, G;

            Doo(int w_, double F_, double G_)
                : w(w_), F(F_), G(G_)
            {
            }
        };

        std::deque<Doo> FG;

        struct SimNode
        {
            int c;
            int64_t sr1, sr2;
            uint64_t shiftFG;

            SimNode(int c_, int64_t sr1_, int64_t sr2_, uint64_t shiftFG_)
                : c(c_), sr1(sr1_), sr2(sr2_), shiftFG(shiftFG_)
            {
            }
        };

        std::deque<SimNode> simnodes;

        void emplace_back(int c_, int64_t sr1_, int64_t sr2_)
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
            int w;
            double E;
        };

        std::unique_ptr<Do[]> E0, E;
        std::unique_ptr<double[]> Ethres;
    };

    const static int ww = 3;

    std::vector<uint8_t> O;
    std::string Oname;
    std::ofstream fout;
    int64_t max_id = 0;
    int Omax;

    CrossGlobalData crossglobaldata;

    MonoDeque<Dot *> dot_sources;

    void apply_memory()
    {
        trn = 0;
        tnn = 1 + nodes.size();
        tsn = targets.size();
        for (auto &node : nodes)
        {
            trn += node.get_trn();
            tnn += node.get_tnn(Omax);
            tsn += node.get_tsn(Omax);
        }
        for (auto &edge : local_crosses)
        {
            trn += edge.get_trn();
            tnn += edge.get_tnn(Omax);
            tsn += edge.get_tsn(Omax);
        }
        for (auto &edge : local_circuits)
        {
            trn += edge.get_trn();
            tnn += edge.get_tnn(Omax);
            tsn += edge.get_tsn(Omax);
        }
        for (auto &edge : global_circuits)
        {
            trn += edge.get_trn();
            tnn += edge.get_tnn(Omax);
            tsn += edge.get_tsn(Omax);
        }
        visits.reset(new bool[tnn]);
        for (int i = 0; i < tnn; ++i)
            visits[i] = false;
        vals.reset(new double[tnn]);
        sources.reset(new double *[tsn]);
        sourcess.reset(new double **[tnn]);
        s_szs.reset(new int[tnn]);
        ids.reset(new int64_t[tnn]);
        ss.reset(new int[trn]);
        ns.reset(new int[trn]);

        bool *fpvisit = visits.get(), *rpvisit = fpvisit + tnn;
        pQvisit = --rpvisit;
        double *fpval = vals.get(), *rpval = fpval + tnn;
        pQval = --rpval;
        double **fpsource = sources.get();
        double ***fpsources = sourcess.get(), ***rpsources = fpsources + tnn;
        pQsources = --rpsources;
        *pQsources = fpsource;
        fpsource += targets.size();
        int *fps_sz = s_szs.get(), *rps_sz = fps_sz + tnn;
        pQs_sz = --rps_sz;
        int64_t *fpid = ids.get(), *rpid = fpid + tnn;
        pQid = --rpid;
        int *fps = ss.get();
        int *fpn = ns.get();
        for (auto &node : nodes)
            node.apply_memory(Omax, fpvisit, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn, rpvisit, rpval, rpsources, rps_sz, rpid);
        for (auto &edge : local_crosses)
            edge.apply_memory(Omax, fpvisit, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn);
        for (auto &edge : local_circuits)
            edge.apply_memory(Omax, fpvisit, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn);
        for (auto &edge : global_circuits)
            edge.apply_memory(Omax, fpvisit, fpval, fpsource, fpsources, fps_sz, fpid, fps, fpn);
    }

    struct SWN
    {
        int64_t s;
        int w;
        int n;
    };

    SWN get_swn(int64_t idx)
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

    Align(boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2short, std::map<std::string, RankVec> &file2rankvec, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2long, std::string mgfile, int Omax_)
        : Graph(vm, file2short, file2rankvec, file2long), fout(mgfile, std::ifstream::binary), Omax(Omax_)
    {
        apply_memory();

        crossglobaldata.E0.reset(new CrossGlobalData::Do[Omax + 1]);
        crossglobaldata.E.reset(new CrossGlobalData::Do[Omax + 1]);
        crossglobaldata.Ethres.reset(new double[Omax + 1]);
    }

    void Mix()
    {
        for (uint64_t i = 0; i < nodes.size(); ++i)
            for (uint64_t w = 0; w <= O.size(); ++w)
                nodes[i].Avals[0][w] = -inf;
        for (auto &scc : sccs)
        {
            for (uint64_t w = 0; w <= O.size(); ++w)
            {
                for (auto node : scc.nodes)
                {
                    if (w == 0)
                        source_max(&node->Bvals[w], {}, {});
                    else
                        source_max(&node->Bvals[w], {&node->Bvals[w - 1], &node->Avals[node->scc_sz - 1][w - 1]}, {node->Bvals[w - 1] + node->ue, node->Avals[node->scc_sz - 1][w - 1] + node->ve});

                    node->updateA0(w, node->Bvals[w], &node->Bvals[w]);
                    if (w == 0)
                        node->updateA0(w, node->pAbarval[0], node->pAbarval);
                }
                for (auto edge : scc.local_circuits)
                    CircuitIteration(w, edge);
                CircuitIterationGlobal(w, scc);
                for (uint64_t l = 1; l < scc.nodes.size(); ++l)
                {
                    for (auto edge : scc.global_circuits)
                        source_max(&edge->D0vals[l - 1][w], {&edge->tail->Avals[l - 1][w]}, {edge->tail->Avals[l - 1][w] + edge->T});
                    for (auto edge : scc.local_circuits)
                    {
                        source_max(&edge->D0vals[l - 1][w], {&edge->tail->Avals[l - 1][w]}, {edge->tail->Avals[l - 1][w] + edge->T});
                        source_max(&edge->DXvals[l - 1][w], {&edge->tail->Avals[l - 1][w], &edge->D0vals[l - 1][w]}, {edge->tail->Avals[l - 1][w] + edge->gfpT[edge->ref_sz], edge->D0vals[l - 1][w] + edge->gf[edge->ref_sz]});
                    }
                    for (auto node : scc.nodes)
                    {
                        node->Avals[l][w] = node->Avals[0][w];
                        for (auto edge : node->in_local_circuits)
                            node->Avals[l][w] = std::max({node->Avals[l][w], edge->D0vals[l - 1][w] + edge->gfm, edge->DXvals[l - 1][w]});
                        for (auto edge : node->in_global_circuits)
                            node->Avals[l][w] = std::max(node->Avals[l][w], edge->D0vals[l - 1][w]);
                        int64_t idx = &node->Avals[l][w] - vals.get();
                        s_szs[idx] = 0;
                        if (node->Avals[0][w] == node->Avals[l][w])
                            sourcess[idx][s_szs[idx]++] = &node->Avals[0][w];
                        for (auto edge : node->in_local_circuits)
                        {
                            if (edge->D0vals[l - 1][w] + edge->gfm == node->Avals[l][w])
                                sourcess[idx][s_szs[idx]++] = &edge->D0vals[l - 1][w];
                            if (edge->DXvals[l - 1][w] == node->Avals[l][w])
                                sourcess[idx][s_szs[idx]++] = &edge->DXvals[l - 1][w];
                        }
                        for (auto edge : node->in_global_circuits)
                            if (edge->D0vals[l - 1][w] == node->Avals[l][w])
                                sourcess[idx][s_szs[idx]++] = &edge->D0vals[l - 1][w];
                    }
                }
            }
            for (auto edge : scc.local_crosses)
                CrossIteration(edge);
            for (auto edge : scc.global_crosses)
                CrossIterationGlobal(edge);
        }
        *pQval = -inf;
        for (auto target : targets)
            *pQval = std::max(*pQval, target->Avals[target->scc_sz - 1][O.size()]);
        s_szs[tnn - 1] = 0;
        for (auto target : targets)
            if (target->Avals[target->scc_sz - 1][O.size()] == *pQval)
                sourcess[tnn - 1][s_szs[tnn - 1]++] = &target->Avals[target->scc_sz - 1][O.size()];
    }

    void CrossIteration(EdgeLocalCross *edge)
    {
        double *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        double **Evals = edge->Evals.get(), **Fvals = edge->Fvals.get(), **Gvals = edge->Gvals.get();
        for (uint64_t s = 0; s <= edge->ref_sz; ++s)
            for (uint64_t w = 0; w <= O.size(); ++w)
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
                    source_max(&Gvals[s][w], {&Fvals[s][w], &Evals[s][w], &Gvals[s - 1][w - 1], &Avals[w]}, {Fvals[s][w], Evals[s][w], Gvals[s - 1][w - 1] + edge->gamma[edge->ref[s - 1]][O[w - 1]], Avals[w] + edge->gfpT[s]});
                edge->head->updateA0(w, Gvals[s][w], &Gvals[s][w]);
            }
    }

    void CrossIterationGlobal(EdgeGlobalCross *edge)
    {
        double *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        std::deque<CrossGlobalData::SimNode> &simnodes = crossglobaldata.simnodes;
        std::unique_ptr<CrossGlobalData::Do[]> &E0 = crossglobaldata.E0, &E = crossglobaldata.E;
        std::unique_ptr<double[]> &Ethres = crossglobaldata.Ethres;
        std::deque<CrossGlobalData::Doo> &FG = crossglobaldata.FG;

        int c = -1;
        int64_t sr1 = 0, sr2 = edge->prankvec->bwt_sz - 1;
        crossglobaldata.emplace_back(c, sr1, sr2);
        for (uint64_t w = 0; w <= O.size(); ++w)
        {
            double tgw = (O.size() > w ? (O.size() - w - 1) * edge->tail->ue + edge->tail->ve : 0);
            E0[w].w = w;
            if (w > 0)
                E0[w].E = std::max(E0[w - 1].E + edge->ue, FG[w - 1].G + edge->ve);
            else
                E0[w].E = -inf;
            double tmp = std::max(E0[w].E, Avals[w] + edge->T);
            if (E0[w].E + edge->ue < tmp + edge->tail->ve && E0[w].E + (O.size() - w) * edge->ue < tmp + tgw)
                E0[w].E = -inf;
            FG.emplace_back(w, -inf, std::max(E0[w].E, Avals[w] + edge->T));
            Ethres[w] = std::max(E0[w].E, FG[w].G + std::min({0.0, edge->tail->ve - edge->ue, tgw - (O.size() - w) * edge->ue}));
            edge->head->updateA0(w, FG[w].G, edge, simnodes[0].sr1, 0);
        }

        while (true)
        {
            do
            {
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
                    sr1 = simnodes.back().sr1;
                    sr2 = simnodes.back().sr2;
                    edge->prankvec->PreRange(sr1, sr2, c);
                } while (sr1 > sr2 && c < 6);
                if (sr1 <= sr2)
                {
                    crossglobaldata.emplace_back(c, sr1, sr2);
                    break;
                }
                else if (FG.size() > simnodes.back().shiftFG)
                    FG.erase(FG.begin() + simnodes.back().shiftFG, FG.end());
            } while (true);
            if (simnodes.empty())
                break;

            uint64_t idx = simnodes[simnodes.size() - 2].shiftFG;
            for (uint64_t e = 0, w = FG[idx].w;; ++e)
            {
                E[e].w = w;
                E[e].E = std::max((e > 0 && E[e - 1].w == w - 1) ? E[e - 1].E + edge->ue : -inf, (simnodes.back().shiftFG < FG.size() && FG.back().w == w - 1) ? FG.back().G + edge->ve : -inf);
                if (E[e].E < Ethres[w])
                    E[e].E = -inf;

                double F_val = idx < simnodes.back().shiftFG && FG[idx].w == w ? std::max(FG[idx].F + edge->uf, FG[idx].G + edge->vf) : -inf;
                double G_val = idx > simnodes[simnodes.size() - 2].shiftFG && FG[idx - 1].w == w - 1 ? FG[idx - 1].G + edge->gamma[simnodes.back().c][O[w - 1]] : -inf;
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

    void CircuitIteration(int64_t w, EdgeLocalCircuit *edge)
    {
        double *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        double *Dvals = edge->Dvals;
        double **Evals = edge->Evals.get(), **F0vals = edge->F0vals.get(), **G0vals = edge->G0vals.get(), **Gvals = edge->Gvals.get();

        if (w > 0)
            source_max(&Dvals[w - 1], {&Avals[w - 1]}, {Avals[w - 1] + edge->T});
        for (uint64_t s = 0; s <= edge->ref_sz; ++s)
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
                source_max(&G0vals[s][w], {&Evals[s][w], &F0vals[s][w], &Gvals[s - 1][w - 1]}, {Evals[s][w], F0vals[s][w], Gvals[s - 1][w - 1] + edge->gamma[edge->ref[s - 1]][O[w - 1]]});
            else
                source_max(&G0vals[s][w], {&Evals[s][w], &F0vals[s][w]}, {Evals[s][w], F0vals[s][w]});
            edge->head->updateA0(w, G0vals[s][w], &G0vals[s][w]);
        }
    }

    void CircuitIterationGlobal(int64_t w, SCC &scc)
    {
        if (w == 0)
            return;

        std::vector<EdgeGlobalCircuit *> edges = scc.global_circuits;
        std::vector<double> tgws(edges.size());
        std::vector<SNC *> jumps(edges.size());

        for (int i = 0; i < edges.size(); ++i)
        {
            if (w == 1)
            {
                edges[i]->sncs.emplace_back(-inf, -inf, -inf, -inf, -1, 5, 0, 0, edges[i]->prankvec->bwt_sz - 1, 1, (SNC *)NULL, (SNC *)NULL);
                for (int c = 2; c <= 6; ++c)
                {
                    int64_t sr1 = edges[i]->sncs.front().sr1;
                    int64_t sr2 = edges[i]->sncs.front().sr2;
                    edges[i]->prankvec->PreRange(sr1, sr2, c);
                    if (sr1 <= sr2)
                    {
                        edges[i]->sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, 1, sr1, sr2, 0, &(edges[i]->sncs.front()), edges[i]->sncs.front().jump);
                        edges[i]->sncs.front().jump = &(edges[i]->sncs.back());
                    }
                }
            }
            else
            {
                for (auto prejump = &(edges[i]->sncs.front()), jump = prejump->jump; jump; jump = prejump->jump)
                {
                    jump->hatG = jump->G0;
                    if (jump->itp != &(edges[i]->sncs.front()) && jump->itp->hatG == -inf && jump->hatG == -inf && jump->E == -inf)
                        prejump->jump = jump->jump;
                    else
                        prejump = jump;
                }
            }

            edges[i]->sncs.front().hatG = std::max(edges[i]->sncs.front().G0, edges[i]->tail->Avals[edges[i]->tail->scc_sz - 1][w - 1] + edges[i]->T);
            edges[i]->sncs.front().E = std::max(edges[i]->sncs.front().hatG + edges[i]->ve, edges[i]->sncs.front().E + edges[i]->ue);
            tgws[i] = (O.size() > w ? (O.size() - w - 1) * edges[i]->tail->ue + edges[i]->tail->ve : 0);
            edges[i]->sncs.front().F0 = -inf;
            edges[i]->sncs.front().G0 = edges[i]->sncs.front().E;
            edges[i]->head->updateA0(w, edges[i]->sncs.front().G0, edges[i], edges[i]->sncs.front().sr1, edges[i]->sncs.front().lambda);

            jumps[i] = edges[i]->sncs.front().jump;
        }

        while (edges.size() > 0)
            for (int i = 0; i < edges.size();)
            {
                double Gthres = std::max(edges[i]->sncs.front().E, edges[i]->tail->Avals[0][w] + edges[i]->T);
                double Ethres = std::max(edges[i]->sncs.front().E, Gthres + std::min(0.0, std::min(edges[i]->tail->ve - edges[i]->ue, tgws[i] - (O.size() - w) * edges[i]->ue)));

                if (jumps[i])
                {
                    jumps[i]->E = std::max(jumps[i]->hatG + edges[i]->ve, jumps[i]->E + edges[i]->ue);
                    if (jumps[i]->E < Ethres)
                        jumps[i]->E = -inf;
                    jumps[i]->F0 = std::max(jumps[i]->itp->G0 + edges[i]->vf, jumps[i]->itp->F0 + edges[i]->uf);
                    jumps[i]->G0 = std::max(std::max(jumps[i]->E, jumps[i]->F0), jumps[i]->itp->hatG + edges[i]->gamma[jumps[i]->c][O[w - 1]]);
                    if (jumps[i]->G0 < Gthres)
                    {
                        jumps[i]->F0 = -inf;
                        jumps[i]->G0 = -inf;
                    }
                    else
                        edges[i]->head->updateA0(w, jumps[i]->G0, edges[i], jumps[i]->sr1, jumps[i]->lambda);
                    if (jumps[i]->G0 != -inf && jumps[i]->hatG == -inf)
                    {
                        if (jumps[i]->cid == 0)
                        {
                            jumps[i]->cid = edges[i]->sncs.size();
                            for (int c = 2; c <= 6; ++c)
                            {
                                int64_t sr1 = jumps[i]->sr1;
                                int64_t sr2 = jumps[i]->sr2;
                                edges[i]->prankvec->PreRange(sr1, sr2, c);
                                if (sr1 <= sr2)
                                {
                                    edges[i]->sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, jumps[i]->lambda + 1, sr1, sr2, 0, jumps[i], (SNC *)NULL);
                                    ++jumps[i]->idl;
                                }
                            }
                        }
                        for (int idi = 0; idi < jumps[i]->idl; ++idi)
                            if (edges[i]->sncs[jumps[i]->cid + idi].E == -inf && edges[i]->sncs[jumps[i]->cid + idi].G0 == -inf)
                            {
                                edges[i]->sncs[jumps[i]->cid + idi].jump = jumps[i]->jump;
                                jumps[i]->jump = &(edges[i]->sncs[jumps[i]->cid + idi]);
                            }
                    }
                    jumps[i] = jumps[i]->jump;
                    
                    ++ i;
                }
                else
                {
                    if (w == O.size())
                        edges[i]->sncs.clear();
                    edges.erase(edges.begin() + i);
                    tgws.erase(tgws.begin() + i);
                    jumps.erase(jumps.begin() + i);
                }
            }
    }

    void outputdot(int64_t M, SWN &swn)
    {
        fout.write((char *)&swn.n, sizeof(Dot::n));
        fout.write((char *)&swn.s, sizeof(Dot::s));
        fout.write((char *)&swn.w, sizeof(Dot::w));
        fout.write((char *)&vals[M], sizeof(Dot::val));
        fout.write((char *)&s_szs[M], sizeof(Dot::s_sz));
        fout.write((char *)&s_szs[M], sizeof(Dot::lambda)); // strange
    }

    void outputdot(Dot *M)
    {
        fout.write((char *)&M->n, sizeof(Dot::n));
        fout.write((char *)&M->s, sizeof(Dot::s));
        fout.write((char *)&M->w, sizeof(Dot::w));
        fout.write((char *)&M->val, sizeof(Dot::val));
        fout.write((char *)&M->s_sz, sizeof(Dot::s_sz));
        fout.write((char *)&M->lambda, sizeof(Dot::lambda));
    }

    void arrangesource(double *source, std::queue<int64_t> &localqueue, int64_t &id)
    {
        int64_t idx = source - vals.get();
        if (!visits[idx])
        {
            ids[idx] = ++id;
            visits[idx] = true;
            localqueue.push(idx);
        }
        fout.write((char *)&ids[idx], sizeof(Dot::id));
    }

    void arrangesource(Dot *source, std::queue<Dot *> &globalqueue, int64_t &id)
    {
        if (!source->visit)
        {
            source->id = ++id;
            source->visit = true;
            globalqueue.push(source);
        }
        fout.write((char *)&source->id, sizeof(Dot::id));
    }

    void GetMinimalGraph()
    {
        int Onsz = Oname.size();
        fout.write((char *)&Onsz, sizeof(Onsz));
        fout.write(Oname.data(), sizeof(char) * Onsz);
        std::queue<Dot *> globalqueue;
        std::queue<int64_t> localqueue;
        std::queue<bool *> visitqueue;
        int64_t id = 0;
        *pQid = id;
        *pQvisit = true;
        localqueue.push(pQval - vals.get());
        while (!globalqueue.empty() || !localqueue.empty())
        {
            if (globalqueue.empty() || !localqueue.empty() && ids[localqueue.front()] < globalqueue.front()->id)
            {
                int64_t M = localqueue.front();
                localqueue.pop();
                visitqueue.push(&visits[M]);
                SWN swn = get_swn(M);
                if (swn.n < -3 && swn.s == 0)
                    GlobalTrack(M, swn, globalqueue, localqueue, id);
                else
                {
                    outputdot(M, swn);
                    for (int i = 0; i < s_szs[M]; ++i)
                        arrangesource(sourcess[M][i], localqueue, id);
                }
            }
            else
            {
                Dot *M = globalqueue.front();
                globalqueue.pop();
                visitqueue.push(&M->visit);
                outputdot(M);
                for (int i = 0; i < M->s_sz; ++i)
                {
                    if (dot_sources[M->fs + i])
                        arrangesource(dot_sources[M->fs + i], globalqueue, id);
                    else
                    {
                        Node *tail = edges[M->n]->tail;
                        arrangesource(&tail->Avals[tail->scc_sz - 1][M->w], localqueue, id);
                    }
                }
            }
        }
        if (id > max_id)
            max_id = id;

        while (!visitqueue.empty())
        {
            *(visitqueue.front()) = false;
            visitqueue.pop();
        }
        for (auto &edge : global_crosses)
            edge.tracktree.clear();
        for (auto &edge : global_circuits)
            edge.tracktree.clear();
        for (auto &node : nodes)
            node.clearAdelta(O.size());
        dot_sources.clear();
    }

    void GlobalTrack(int64_t M, SWN &swn, std::queue<Dot *> &globalqueue, std::queue<int64_t> &localqueue, int64_t &id)
    {
        auto &node = nodes[Dot::nidx_trans(swn.n)];
        s_szs[M] = node.AdeltaGlobal[swn.w].size() + node.AdeltaDot[swn.w].size();
        outputdot(M, swn);
        for (double *source : node.AdeltaDot[swn.w])
            arrangesource(source, localqueue, id);

        for (auto &globalsuffix : node.AdeltaGlobal[swn.w])
        {
            EdgeGlobal *edge = globalsuffix.edge;
            double *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
            TrackTree &tracktree = edge->tracktree;
            std::deque<TrackTree::TrackNode> &tracknodes = tracktree.tracknodes;
            std::deque<Dot> &dots = tracktree.dots;
            uint8_t *longref = file2long[edge->name].first.get();
            std::ifstream &SAfin = file2SA[edge->name];
            int64_t start;
            SAfin.seekg(globalsuffix.start * sizeof(uint64_t));
            SAfin.read((char*)&start, sizeof(uint64_t));
            if (tracktree.tracknodes.empty())
            {
                tracktree.W = O.size();
                tracktree.emplace_back(-1, edge->n, edge->prankvec->bwt_sz - 1 - start - globalsuffix.lambda, 0);
            }
            tracktree.setidx(0);
            if (edge->head->scc_id != edge->tail->scc_id)
            {
                for (int w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w == 0)
                        source_max(dots[tracktree.shiftE + w], {}, {});
                    else
                        source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                    source_max(dots[tracktree.shiftF + w], {}, {});
                    source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], NULL}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, Avals[w] + edge->T});
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                for (int i = globalsuffix.lambda - 1; i >= 0; --i)
                {
                    uint8_t c = longref[start + i];
                    if (tracknodes[tracktree.idx].cidxs[c - 2] < 0)
                    {
                        tracktree.emplace_back(-1, edge->n, edge->prankvec->bwt_sz - 1 - start - i, dots[tracktree.shiftE].lambda + 1);
                        tracknodes[tracktree.idx].cidxs[c - 2] = tracknodes.size() - 1;
                    }
                    tracktree.setidx(tracknodes[tracktree.idx].cidxs[c - 2]);
                    for (int w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
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
                for (int w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
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
                for (int i = globalsuffix.lambda - 1; i >= 0; --i)
                {
                    uint8_t c = longref[start + i];
                    if (tracknodes[tracktree.idx].cidxs[c - 2] < 0)
                    {
                        tracktree.emplace_back(-1, edge->n, edge->prankvec->bwt_sz - 1 - start - i, dots[tracktree.shiftE].lambda + 1);
                        tracknodes[tracktree.idx].cidxs[c - 2] = tracknodes.size() - 1;
                    }
                    tracktree.setidx(tracknodes[tracktree.idx].cidxs[c - 2]);
                    for (int w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
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
            arrangesource(&dots[tracktree.shiftG + swn.w], globalqueue, id);
        }
    }

    void source_max(double *val, std::initializer_list<double *> sources, std::initializer_list<double> source_vals)
    {
        int64_t idx = val - vals.get();
        s_szs[idx] = 0;
        if (source_vals.begin() == source_vals.end())
            *val = -inf;
        else
        {
            *val = std::max(source_vals);
            auto sources_it = sources.begin();
            for (auto source_vals_it = source_vals.begin(); source_vals_it != source_vals.end(); ++source_vals_it, ++sources_it)
                if (*source_vals_it == *val)
                    sourcess[idx][s_szs[idx]++] = *sources_it;
        }
    }

    void source_max(Dot &dot, std::initializer_list<Dot *> adrs, std::initializer_list<double> vals)
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
