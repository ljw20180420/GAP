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

struct SNC
{
    double E;
    double F0;
    double G0;
    double hatG;
    int8_t c;
    int8_t idl;
    int lambda;
    size_t sr1;
    size_t sr2;
    size_t cid; // 0 means not try to add child yet
    SNC *itp;
    SNC *jump;

    SNC(double E_, double F0_, double G0_, double hatG_, int8_t c_, int8_t idl_, int lambda_, size_t sr1_, size_t sr2_, size_t cid_, SNC *itp_, SNC *jump_)
        : E(E_), F0(F0_), G0(G0_), hatG(hatG_), c(c_), idl(idl_), lambda(lambda_), sr1(sr1_), sr2(sr2_), cid(cid_), itp(itp_), jump(jump_)
    {
    }
};

struct Align : Memory, Graph
{

    struct CrossGlobalData
    {
        struct SimNode
        {
            struct Doo
            {
                int w;
                double F;
                double G;

                Doo(int w_, double F_, double G_)
                    : w(w_), F(F_), G(G_)
                {
                }
            };
            int c;
            size_t sr1;
            size_t sr2;
            std::deque<Doo> FG;

            SimNode(int c_, size_t sr1_, size_t sr2_)
                : c(c_), sr1(sr1_), sr2(sr2_)
            {
            }

            void reset(int c_, size_t sr1_, size_t sr2_)
            {
                c = c_;
                sr1 = sr1_;
                sr2 = sr2_;
            }
        };

        struct Do
        {
            int w;
            double E;

            Do(int w_, double E_)
            {
                w = w_;
                E = E_;
            }
        };

        std::deque<Do> E0, E;
        std::deque<SimNode> simnodes;
        std::deque<std::deque<std::pair<size_t, int>> *> Cdelta;
        std::deque<double> Ethres;
    };

    const static int ww = 3;

    std::string O;
    std::string Oname;
    Dot Q;
    std::vector<int> target_scc_szes;
    std::vector<int> node_scc_szes;
    std::ofstream fout;
    size_t max_id = 0;

    CrossGlobalData crossglobaldata;

    void alloc_initial(Dot **&ptr, int n, int S)
    {
        ptr = heap_alloc<Dot *>(S + 1);
        for (int s = 0; s <= S; ++s)
            alloc_initial(ptr[s], n, s);
    }

    void alloc_initial(Dot *&ptr, int n, int s)
    {
        ptr = heap_alloc<Dot>(O.size() + 1);
        for (size_t w = 0; w <= O.size(); ++w)
        {
            ptr[w].n = n;
            ptr[w].s = s;
            ptr[w].w = w;
        }
    }

    Align(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel, size_t chunk_sz_, std::string file) : Graph(argc, argv, file2seq, file2browheel), fout(file, std::ifstream::binary)
    {
        Initial(chunk_sz_);
        for (auto &node : nodes)
        {
            node.Abar.n = -2;
            node.Abar.s = 0;
            node.Abar.w = 0;
            node.Abar.val = -inf;
            node.Abar.s_sz = 0;
            node.Abar.sources = NULL;
        }
        for (auto root : roots)
            root->Abar.val = 0;
        node_scc_szes.resize(nodes.size());
        target_scc_szes.resize(targets.size());
        for (auto &scc : sccs)
            for (auto node : scc.nodes)
            {
                for (size_t i = 0; i < targets.size(); ++i)
                    if (targets[i] == node)
                    {
                        target_scc_szes[i] = scc.nodes.size();
                        break;
                    }
                for (size_t i = 0; i < nodes.size(); ++i)
                    if (&nodes[i] == node)
                        node_scc_szes[i] = scc.nodes.size();
            }
        Q.n = -3;
        Q.s = 0;
        Q.w = 0;
    }

    void Mix()
    {
        if (O.size() + 1 > crossglobaldata.Cdelta.size())
        {
            int rsz = std::max(O.size() + 1, 2 * crossglobaldata.Cdelta.size());
            crossglobaldata.Cdelta.resize(rsz);
            crossglobaldata.Ethres.resize(rsz);
            for (auto &node : nodes)
            {
                node.AdeltaCross.resize(rsz);
                node.AdeltaCircuit.resize(rsz);
            }
        }
        for (auto &scc : sccs)
        {
            for (auto node : scc.nodes)
            {
                alloc_initial(node->A, nodes.data() - node - 4, scc.nodes.size() - 1);
                alloc_initial(node->B, -1, 0);
            }
            for (auto edge : scc.local_circuits)
            {
                alloc_initial(edge->D, edge->n, 0);
                alloc_initial(edge->E, edge->n, edge->nameseq.seq.size());
                alloc_initial(edge->F0, edge->n, edge->nameseq.seq.size());
                alloc_initial(edge->G0, edge->n, edge->nameseq.seq.size());
                alloc_initial(edge->G, edge->n, edge->nameseq.seq.size());
                alloc_initial(edge->D0, edge->n, scc.nodes.size() - 1);
                alloc_initial(edge->DX, edge->n, scc.nodes.size() - 1);
            }
            for (auto edge : scc.global_circuits)
            {
                alloc_initial(edge->D0, edge->n, scc.nodes.size() - 1);
                edge->C.resize(O.size() + 1);
            }
            for (size_t w = 0; w <= O.size(); ++w)
            {
                for (auto node : scc.nodes)
                    if (w == 0)
                        source_max(node->B[w], {}, {});
                    else
                        source_max(node->B[w], {&node->B[w - 1], &node->A[scc.nodes.size() - 1][w - 1]}, {node->B[w - 1].val + node->ue, node->A[scc.nodes.size() - 1][w - 1].val + node->ve});
                for (auto edge : scc.local_circuits)
                    CircuitIteration(w, edge, scc.nodes.size());
                for (auto edge : scc.global_circuits)
                    CircuitIterationGlobal(w, edge, scc.nodes.size());
                for (auto node : scc.nodes)
                {
                    std::vector<Dot *> adrs{&node->B[w]};
                    std::vector<double> vals{node->B[w].val};
                    for (auto edge : node->in_local_crosses)
                        for (size_t s = 0; s <= edge->nameseq.seq.size(); ++s)
                        {
                            adrs.push_back(&edge->G[s][w]);
                            vals.push_back(edge->G[s][w].val + (s == edge->nameseq.seq.size() ? 0 : edge->vfm + (edge->nameseq.seq.size() - s - 1) * edge->ufm));
                        }
                    for (auto edge : node->in_local_circuits)
                        for (size_t s = 0; s <= edge->nameseq.seq.size(); ++s)
                        {
                            adrs.push_back(&edge->G0[s][w]);
                            vals.push_back(edge->G0[s][w].val + (s == edge->nameseq.seq.size() ? 0 : edge->vfm + (edge->nameseq.seq.size() - s - 1) * edge->ufm));
                        }
                    if (w == 0)
                    {
                        adrs.push_back(&node->Abar);
                        vals.push_back(node->Abar.val);
                    }
                    source_max(node->A[0][w], adrs.begin(), vals.begin(), vals.end());
                    double local_max = node->A[0][w].val;
                    for (auto edge : node->in_global_crosses)
                        if (edge->C[w] > node->A[0][w].val)
                            node->A[0][w].val = edge->C[w];
                    for (auto edge : node->in_global_circuits)
                        if (w > 0 && edge->C[w] > node->A[0][w].val)
                            node->A[0][w].val = edge->C[w];
                    if (local_max < node->A[0][w].val)
                        node->A[0][w].s_sz = 0;
                    for (auto edge : node->in_global_crosses)
                        if (edge->C[w] < node->A[0][w].val)
                            node->AdeltaCross[w].erase(edge);
                    for (auto edge : node->in_global_circuits)
                        if (w > 0 && edge->C[w] < node->A[0][w].val)
                            node->AdeltaCircuit[w].erase(edge);
                }
                for (size_t l = 1; l < scc.nodes.size(); ++l)
                {
                    for (auto edge : scc.global_circuits)
                        source_max(edge->D0[l][w], {&edge->tail->A[l - 1][w]}, {edge->tail->A[l - 1][w].val + edge->T});
                    for (auto edge : scc.local_circuits)
                    {
                        source_max(edge->D0[l][w], {&edge->tail->A[l - 1][w]}, {edge->tail->A[l - 1][w].val + edge->T});
                        source_max(edge->DX[l][w], {&edge->tail->A[l - 1][w], &edge->D0[l][w]}, {edge->tail->A[l - 1][w].val + edge->vfp + (edge->nameseq.seq.size() - 1) * edge->ufp + edge->T, edge->D0[l][w].val + edge->vf + (edge->nameseq.seq.size() - 1) * edge->uf});
                    }
                    for (auto node : scc.nodes)
                    {
                        std::vector<Dot *> adrs{&node->A[0][w]};
                        std::vector<double> vals{node->A[0][w].val};
                        for (auto edge : node->in_local_circuits)
                        {
                            adrs.push_back(&edge->D0[l][w]);
                            vals.push_back(edge->D0[l][w].val + edge->vfm + (edge->nameseq.seq.size() - 1) * edge->ufm);
                            adrs.push_back(&edge->DX[l][w]);
                            vals.push_back(edge->DX[l][w].val);
                        }
                        for (auto edge : node->in_global_circuits)
                        {
                            adrs.push_back(&edge->D0[l][w]);
                            vals.push_back(edge->D0[l][w].val);
                        }
                        source_max(node->A[l][w], adrs.begin(), vals.begin(), vals.end());
                    }
                }
            }
            for (auto edge : scc.local_crosses)
            {
                alloc_initial(edge->E, edge->n, edge->nameseq.seq.size());
                alloc_initial(edge->F, edge->n, edge->nameseq.seq.size());
                alloc_initial(edge->G, edge->n, edge->nameseq.seq.size());
                CrossIteration(edge, scc.nodes.size());
            }
            for (auto edge : scc.global_crosses)
            {
                edge->C.resize(O.size() + 1);
                CrossIterationGlobal(edge, scc.nodes.size());
            }
        }
        std::vector<Dot *> adrs;
        std::vector<double> vals;
        for (size_t n = 0; n < targets.size(); ++n)
        {
            adrs.push_back(&targets[n]->A[target_scc_szes[n] - 1][O.size()]);
            vals.push_back(targets[n]->A[target_scc_szes[n] - 1][O.size()].val);
        }
        source_max(Q, adrs.begin(), vals.begin(), vals.end());
    }

    void CrossIteration(EdgeLocalCross *edge, int scc_sz)
    {
        Dot **A = edge->tail->A;
        Dot **E = edge->E;
        Dot **F = edge->F;
        Dot **G = edge->G;
        for (size_t s = 0; s <= edge->nameseq.seq.size(); ++s)
            for (size_t w = 0; w <= O.size(); ++w)
            {
                if (w == 0)
                    source_max(E[s][w], {}, {});
                else
                    source_max(E[s][w], {&E[s][w - 1], &G[s][w - 1]}, {E[s][w - 1].val + edge->ue, G[s][w - 1].val + edge->ve});
                if (s == 0)
                    source_max(F[s][w], {}, {});
                else
                    source_max(F[s][w], {&F[s - 1][w], &G[s - 1][w]}, {F[s - 1][w].val + edge->uf, G[s - 1][w].val + edge->vf});
                if (s == 0 || w == 0)
                    source_max(G[s][w], {&F[s][w], &E[s][w], &A[scc_sz - 1][w]}, {F[s][w].val, E[s][w].val, A[scc_sz - 1][w].val + (s == 0 ? 0 : edge->vfp + (s - 1) * edge->ufp) + edge->T});
                else
                    source_max(G[s][w], {&F[s][w], &E[s][w], &G[s - 1][w - 1], &A[scc_sz - 1][w]}, {F[s][w].val, E[s][w].val, G[s - 1][w - 1].val + edge->gamma[BroWheel::base2int(edge->nameseq.seq[s - 1])][BroWheel::base2int(O[w - 1])], A[scc_sz - 1][w].val + (s == 0 ? 0 : edge->vfp + (s - 1) * edge->ufp) + edge->T});
            }
    }

    void CrossIterationGlobal(EdgeGlobalCross *edge, int scc_sz)
    {
        Dot *A = edge->tail->A[scc_sz - 1];
        auto &simnodes = crossglobaldata.simnodes;
        auto &E0 = crossglobaldata.E0;
        auto &E = crossglobaldata.E;
        auto &Cdelta = crossglobaldata.Cdelta;
        auto &Ethres = crossglobaldata.Ethres;
        for (int w = 0; w <= O.size(); ++w)
            Cdelta[w] = &edge->head->AdeltaCross[w][edge];

        int lambda = 0;
        int c = -1;
        int64_t sr1 = 0, sr2 = edge->browheel.sequence.size() - 1;
        if (lambda < simnodes.size())
            simnodes[lambda].reset(c, sr1, sr2);
        else
            simnodes.emplace_back(c, sr1, sr2);
        for (size_t w = 0; w <= O.size(); ++w)
        {
            double tgw = (O.size() > w ? (O.size() - w - 1) * edge->tail->ue + edge->tail->ve : 0);
            E0.emplace_back(w, w > 0 ? std::max(E0[w - 1].E + edge->ue, simnodes[0].FG[w - 1].G + edge->ve) : -inf);
            if (E0[w].E + edge->ue < std::max(E0[w].E, A[w].val + edge->T) + edge->tail->ve && E0[w].E + (O.size() - w) * edge->ue < std::max(E0[w].E, A[w].val + edge->T) + tgw)
                E0[w].E = -inf;
            simnodes[lambda].FG.emplace_back(w, -inf, std::max(E0[w].E, A[w].val + edge->T));
            Ethres[w] = std::max(E0[w].E, simnodes[0].FG[w].G + std::min(0.0, std::min(edge->tail->ve - edge->ue, tgw - (O.size() - w) * edge->ue)));
            edge->C[w] = simnodes[lambda].FG[w].G;
            Cdelta[w]->emplace_back(simnodes[lambda].sr1, lambda);
        }

        size_t totalcount = 0; // debug
        while (true)
        {
            do
            {
                if (!simnodes[lambda].FG.empty())
                    c = 1;
                else
                {
                    do
                    {
                        c = simnodes[lambda].c;
                        simnodes[lambda].FG.clear();
                        --lambda;
                    } while (c == 6 && lambda >= 0);
                    if (lambda < 0)
                        break;
                }
                do
                {
                    ++c;
                    sr1 = simnodes[lambda].sr1;
                    sr2 = simnodes[lambda].sr2;
                    edge->browheel.PreRange(sr1, sr2, c);
                } while (sr1 > sr2 && c < 6);
                if (sr1 <= sr2)
                {
                    ++lambda;
                    if (lambda < simnodes.size())
                        simnodes[lambda].reset(c, sr1, sr2);
                    else
                        simnodes.emplace_back(c, sr1, sr2);
                    break;
                }
                else
                    simnodes[lambda].FG.clear();
            } while (true);
            if (lambda < 0)
                break;

            ++totalcount; // debug
            for (int j = 0, w = simnodes[lambda - 1].FG[j].w;;)
            {
                E.emplace_back(w, std::max((!E.empty() && E.back().w == w - 1) ? E.back().E + edge->ue : -inf, (!simnodes[lambda].FG.empty() && simnodes[lambda].FG.back().w == w - 1) ? simnodes[lambda].FG.back().G + edge->ve : -inf));
                if (E.back().E < Ethres[w])
                    E.back().E = -inf;

                double F_val = j < simnodes[lambda - 1].FG.size() && simnodes[lambda - 1].FG[j].w == w ? std::max(simnodes[lambda - 1].FG[j].F + edge->uf, simnodes[lambda - 1].FG[j].G + edge->vf) : -inf;
                double G_val = j > 0 && simnodes[lambda - 1].FG[j - 1].w == w - 1 ? simnodes[lambda - 1].FG[j - 1].G + edge->gamma[simnodes[lambda].c][BroWheel::base2int(O[w - 1])] : -inf;
                G_val = std::max(std::max(G_val, E.back().E), F_val);
                if (G_val >= simnodes[0].FG[w].G)
                {
                    simnodes[lambda].FG.emplace_back(w, F_val, G_val);
                    if (G_val >= edge->C[w])
                    {
                        if (G_val > edge->C[w])
                        {
                            edge->C[w] = G_val;
                            Cdelta[w]->clear();
                        }
                        Cdelta[w]->emplace_back(simnodes[lambda].sr1, lambda);
                    }
                }
                if (j < simnodes[lambda - 1].FG.size() && simnodes[lambda - 1].FG[j].w <= w)
                    ++j;
                if (w < O.size() && (simnodes[lambda - 1].FG[j - 1].w == w || (!simnodes[lambda].FG.empty() && simnodes[lambda].FG.back().w == w) || (!E.empty() && E.back().w == w && E.back().E > -inf)))
                    ++w;
                else if (j < simnodes[lambda - 1].FG.size())
                    w = simnodes[lambda - 1].FG[j].w;
                else
                    break;
            }

            E.clear();
        }

        E0.clear();

        std::cerr << totalcount << '\n'; // debug
    }

    void CircuitIteration(size_t w, EdgeLocalCircuit *edge, int scc_sz)
    {
        Dot **A = edge->tail->A;
        Dot *D = edge->D;
        Dot **E = edge->E;
        Dot **F0 = edge->F0;
        Dot **G0 = edge->G0;
        Dot **G = edge->G;
        if (w > 0)
            source_max(D[w - 1], {&A[scc_sz - 1][w - 1]}, {A[scc_sz - 1][w - 1].val + edge->T});
        for (size_t s = 0; s <= edge->nameseq.seq.size(); ++s)
        {
            if (w > 0)
            {
                if (s > 0)
                    source_max(G[s][w - 1], {&G0[s][w - 1], &D[w - 1], &A[scc_sz - 1][w - 1]}, {G0[s][w - 1].val, D[w - 1].val + edge->vf + (s - 1) * edge->uf, A[scc_sz - 1][w - 1].val + edge->vfp + (s - 1) * edge->ufp + edge->T});
                else
                    source_max(G[s][w - 1], {&G0[s][w - 1], &D[w - 1]}, {G0[s][w - 1].val, D[w - 1].val});
            }
            if (w == 0)
                source_max(E[s][w], {}, {});
            else
                source_max(E[s][w], {&E[s][w - 1], &G[s][w - 1]}, {E[s][w - 1].val + edge->ue, G[s][w - 1].val + edge->ve});
            if (s == 0)
                source_max(F0[s][w], {}, {});
            else
                source_max(F0[s][w], {&F0[s - 1][w], &G0[s - 1][w]}, {F0[s - 1][w].val + edge->uf, G0[s - 1][w].val + edge->vf});
            if (s > 0 && w > 0)
                source_max(G0[s][w], {&E[s][w], &F0[s][w], &G[s - 1][w - 1]}, {E[s][w].val, F0[s][w].val, G[s - 1][w - 1].val + edge->gamma[BroWheel::base2int(edge->nameseq.seq[s - 1])][BroWheel::base2int(O[w - 1])]});
            else
                source_max(G0[s][w], {&E[s][w], &F0[s][w]}, {E[s][w].val, F0[s][w].val});
        }
    }

    void CircuitIterationGlobal(size_t w, EdgeGlobalCircuit *edge, int scc_sz)
    {
        if (w == 0)
            return;

        Dot *A = edge->tail->A[scc_sz - 1], *B = edge->tail->B;
        std::deque<std::pair<size_t, int>> * Cdelta = &edge->head->AdeltaCircuit[w][edge];
        std::deque<SNC> &sncs = edge->sncs;

        if (w == 1)
        {
            sncs.emplace_back(-inf, -inf, -inf, -inf, -1, 5, 0, 0, edge->browheel.sequence.size() - 1, 1, (SNC *)NULL, (SNC *)NULL);
            for (int c = 2; c <= 6; ++c)
            {
                int64_t sr1 = sncs.front().sr1;
                int64_t sr2 = sncs.front().sr2;
                edge->browheel.PreRange(sr1, sr2, c);
                if (sr1 <= sr2)
                {
                    sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, 1, sr1, sr2, 0, &sncs.front(), sncs.front().jump);
                    sncs.front().jump = &sncs.back();
                }
            }
        }
        else
        {
            for (auto prejump = &sncs.front(), jump = prejump->jump; jump; jump = prejump->jump)
            {
                jump->hatG = jump->G0;
                if (jump->itp != &sncs.front() && jump->itp->hatG == -inf && jump->hatG == -inf && jump->E == -inf)
                    prejump->jump = jump->jump;
                else
                    prejump = jump;
            }
        }

        sncs.front().hatG = std::max(sncs.front().G0, A[w - 1].val + edge->T);
        sncs.front().E = std::max(sncs.front().hatG + edge->ve, sncs.front().E + edge->ue);
        double tgw = (O.size() > w ? (O.size() - w - 1) * edge->tail->ue + edge->tail->ve : 0);
        if (sncs.front().E + edge->ue < std::max(sncs.front().E, B[w].val + edge->T) + edge->tail->ve && sncs.front().E + (O.size() - w) * edge->ue < std::max(sncs.front().E, B[w].val + edge->T) + tgw)
            sncs.front().E = -inf;
        double Gthres = std::max(sncs.front().E, B[w].val + edge->T);
        double Ethres = std::max(sncs.front().E, Gthres + std::min(0.0, std::min(edge->tail->ve - edge->ue, tgw - (O.size() - w) * edge->ue)));
        sncs.front().F0 = -inf;
        sncs.front().G0 = sncs.front().E;
        edge->C[w] = sncs.front().G0;
        Cdelta->emplace_back(sncs.front().sr1, sncs.front().lambda);

        for (auto jump = sncs.front().jump; jump; jump = jump->jump)
        {
            jump->E = std::max(jump->hatG + edge->ve, jump->E + edge->ue);
            if (jump->E < Ethres)
                jump->E = -inf;
            jump->F0 = std::max(jump->itp->G0 + edge->vf, jump->itp->F0 + edge->uf);
            jump->G0 = std::max(std::max(jump->E, jump->F0), jump->itp->hatG + edge->gamma[jump->c][BroWheel::base2int(O[w - 1])]);
            if (jump->G0 < Gthres)
            {
                jump->F0 = -inf;
                jump->G0 = -inf;
            }
            else if (jump->G0 >= edge->C[w])
            {
                if (jump->G0 > edge->C[w])
                {
                    edge->C[w] = jump->G0;
                    Cdelta->clear();
                }
                Cdelta->emplace_back(jump->sr1, jump->lambda);
            }
            if (jump->G0 != -inf && jump->hatG == -inf)
            {
                if (jump->cid == 0)
                {
                    jump->cid = sncs.size();
                    for (int c = 2; c <= 6; ++c)
                    {
                        int64_t sr1 = jump->sr1;
                        int64_t sr2 = jump->sr2;
                        edge->browheel.PreRange(sr1, sr2, c);
                        if (sr1 <= sr2)
                        {
                            sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, jump->lambda + 1, sr1, sr2, 0, jump, (SNC *)NULL);
                            ++jump->idl;
                        }
                    }
                }
                for (int idi = 0; idi < jump->idl; ++idi)
                    if (sncs[jump->cid + idi].E == -inf && sncs[jump->cid + idi].G0 == -inf)
                    {
                        sncs[jump->cid + idi].jump = jump->jump;
                        jump->jump = &sncs[jump->cid + idi];
                    }
            }
        }
        if (w == O.size())
            sncs.clear();
;    }

    void GetMinimalGraph()
    {
        int Onsz = Oname.size();
        fout.write((char *)&Onsz, sizeof(Onsz));
        fout.write(Oname.data(), sizeof(char) * Onsz);
        std::queue<Dot *> dots;
        size_t id = 0;
        Q.id = id;
        Q.s_sz = -Q.s_sz - 1;
        dots.push(&Q);
        while (!dots.empty())
        {
            Dot *M = dots.front();
            dots.pop();
            M->s_sz = -M->s_sz - 1;
            if (M->n < -3 && M->s == 0)
                GlobalTrack(M);
            fout.write((char *)&M->n, sizeof(Dot::n));
            fout.write((char *)&M->s, sizeof(Dot::s));
            fout.write((char *)&M->w, sizeof(Dot::w));
            fout.write((char *)&M->val, sizeof(Dot::val));
            fout.write((char *)&M->s_sz, sizeof(Dot::s_sz));
            fout.write((char *)&M->lambda, sizeof(Dot::lambda));
            for (int i = 0; i < M->s_sz; ++i)
            {
                if (M->sources[i]->s_sz >= 0)
                {
                    M->sources[i]->id = ++id;
                    M->sources[i]->s_sz = -M->sources[i]->s_sz - 1;
                    dots.push(M->sources[i]);
                }
                fout.write((char *)&M->sources[i]->id, sizeof(M->sources[i]->id));
            }
        }
        for (auto &edge : global_crosses)
            edge.tracknodes.clear();
        for (auto &edge : global_circuits)
            edge.tracknodes.clear();
        max_id = id > max_id ? id : max_id;

        for (auto &node : nodes)
            for (int w = 0; w <= O.size(); ++w)
            {
                node.AdeltaCross[w].clear();
                node.AdeltaCircuit[w].clear();
            }
        clear();
    }

    void GlobalTrack(Dot *M)
    {
        std::deque<Dot *> add_sources;
        for (auto &pair : nodes[-M->n - 4].AdeltaCross[M->w])
        {
            Dot **A = pair.first->tail->A;
            auto &tracknodes = pair.first->tracknodes;
            for (auto &suf : pair.second)
            {
                size_t start = pair.first->browheel.SimSuffix(suf.first);
                if (tracknodes.empty())
                    tracknodes.emplace_back(this, (TrackNode *)NULL, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size() - 1 - start - suf.second, O.size(), 0);
                auto tracknode = &tracknodes.front();
                for (int w = tracknode->tau; w <= M->w; ++w)
                {
                    if (w == 0)
                        source_max(tracknode->E[w], {}, {});
                    else
                        source_max(tracknode->E[w], {&tracknode->E[w - 1], &tracknode->G[w - 1]}, {tracknode->E[w - 1].val + pair.first->ue, tracknode->G[w - 1].val + pair.first->ve});
                    source_max(tracknode->F[w], {}, {});
                    source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &A[node_scc_szes[pair.first->tail - nodes.data()] - 1][w]}, {tracknode->E[w].val, tracknode->F[w].val, A[node_scc_szes[pair.first->tail - nodes.data()] - 1][w].val + pair.first->T});
                }
                tracknode->tau = std::max(tracknode->tau, M->w + 1);
                for (int i = suf.second - 1; i >= 0; --i)
                {
                    int c = pair.first->browheel.sequence(start + i);
                    if (!tracknode->itcs[c - 2])
                    {
                        tracknodes.emplace_back(this, tracknode, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size() - 1 - start - i, O.size(), tracknode->lambda + 1);
                        tracknode->itcs[c - 2] = &tracknodes.back();
                    }
                    tracknode = tracknode->itcs[c - 2];
                    for (int w = tracknode->tau; w <= M->w; ++w)
                    {
                        if (w == 0)
                            source_max(tracknode->E[w], {}, {});
                        else
                            source_max(tracknode->E[w], {&tracknode->E[w - 1], &tracknode->G[w - 1]}, {tracknode->E[w - 1].val + pair.first->ue, tracknode->G[w - 1].val + pair.first->ve});
                        source_max(tracknode->F[w], {&tracknode->itp->F[w], &tracknode->itp->G[w]}, {tracknode->itp->F[w].val + pair.first->uf, tracknode->itp->G[w].val + pair.first->vf});
                        if (w == 0)
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w]}, {tracknode->E[w].val, tracknode->F[w].val});
                        else
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &tracknode->itp->G[w - 1]}, {tracknode->E[w].val, tracknode->F[w].val, tracknode->itp->G[w - 1].val + pair.first->gamma[c][BroWheel::base2int(O[w - 1])]});
                    }
                    tracknode->tau = std::max(tracknode->tau, M->w + 1);
                }
                add_sources.push_back(&tracknode->G[M->w]);
            }
        }
        for (auto &pair : nodes[-M->n - 4].AdeltaCircuit[M->w])
        {
            Dot **A = pair.first->tail->A;
            auto &tracknodes = pair.first->tracknodes;
            for (auto &suf : pair.second)
            {
                size_t start = pair.first->browheel.SimSuffix(suf.first);
                if (tracknodes.empty())
                {
                    alloc_initial(pair.first->G, pair.first->n, pair.first->browheel.sequence.size() - 1 - start - suf.second);
                    for (size_t w = 0; w <= O.size(); ++w)
                        pair.first->G[w].lambda = 0;
                    tracknodes.emplace_back(this, (TrackNode *)NULL, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size() - 1 - start - suf.second, O.size(), 0);
                }
                auto tracknode = &tracknodes.front();
                for (int w = tracknode->tau; w <= M->w; ++w)
                {
                    if (w > 0)
                        source_max(pair.first->G[w - 1], {&tracknode->G[w - 1], &A[node_scc_szes[pair.first->tail - nodes.data()] - 1][w - 1]}, {tracknode->G[w - 1].val, A[node_scc_szes[pair.first->tail - nodes.data()] - 1][w - 1].val + pair.first->T});
                    if (w == 0)
                        source_max(tracknode->E[w], {}, {});
                    else
                        source_max(tracknode->E[w], {&tracknode->E[w - 1], &pair.first->G[w - 1]}, {tracknode->E[w - 1].val + pair.first->ue, pair.first->G[w - 1].val + pair.first->ve});
                    source_max(tracknode->F[w], {}, {});
                    source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w]}, {tracknode->E[w].val, tracknode->F[w].val});
                }
                tracknode->tau = std::max(tracknode->tau, M->w + 1);
                for (int i = suf.second - 1; i >= 0; --i)
                {
                    int c = pair.first->browheel.sequence(start + i);
                    if (!tracknode->itcs[c - 2])
                    {
                        tracknodes.emplace_back(this, tracknode, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size() - 1 - start - i, O.size(), tracknode->lambda + 1);
                        tracknode->itcs[c - 2] = &tracknodes.back();
                    }
                    tracknode = tracknode->itcs[c - 2];
                    for (int w = tracknode->tau; w <= M->w; ++w)
                    {
                        if (w == 0)
                            source_max(tracknode->E[w], {}, {});
                        else
                            source_max(tracknode->E[w], {&tracknode->E[w - 1], &tracknode->G[w - 1]}, {tracknode->E[w - 1].val + pair.first->ue, tracknode->G[w - 1].val + pair.first->ve});
                        source_max(tracknode->F[w], {&tracknode->itp->F[w], &tracknode->itp->G[w]}, {tracknode->itp->F[w].val + pair.first->uf, tracknode->itp->G[w].val + pair.first->vf});
                        if (w == 0)
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w]}, {tracknode->E[w].val, tracknode->F[w].val});
                        else if (tracknode->itp == &tracknodes.front())
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &pair.first->G[w - 1]}, {tracknode->E[w].val, tracknode->F[w].val, pair.first->G[w - 1].val + pair.first->gamma[c][BroWheel::base2int(O[w - 1])]});
                        else
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &tracknode->itp->G[w - 1]}, {tracknode->E[w].val, tracknode->F[w].val, tracknode->itp->G[w - 1].val + pair.first->gamma[c][BroWheel::base2int(O[w - 1])]});
                    }
                    tracknode->tau = std::max(tracknode->tau, M->w + 1);
                }
                add_sources.push_back(&tracknode->G[M->w]);
            }
        }
        Dot **old_sources = M->sources;
        int old_s_sz = M->s_sz;
        M->s_sz += add_sources.size();
        M->sources = heap_alloc<Dot *>(M->s_sz);
        auto iter = add_sources.begin();
        for (int i = 0; i < M->s_sz; ++i)
            if (i < old_s_sz)
                M->sources[i] = old_sources[i];
            else
            {
                M->sources[i] = *iter;
                ++iter;
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
        dot.s_sz = 0;
        for (Itv itv = itv_b; itv != itv_e; ++itv)
        {
            if (*itv > dot.val)
                dot.val = *itv, dot.s_sz = 1;
            else if (*itv == dot.val)
                ++dot.s_sz;
        }
        dot.sources = heap_alloc<Dot *>(dot.s_sz);
        int j = 0;
        for (Itv itv = itv_b; itv != itv_e; ++itv, ++ita_b)
            if (*itv == dot.val)
                dot.sources[j++] = *ita_b;
    }
};

#endif
