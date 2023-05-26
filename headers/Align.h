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
        std::deque<double> C, Ethres;
    };

    const static int ww = 3;

    std::string O;
    std::string Oname;
    Dot Q;
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

    void dot_initial(std::vector<std::deque<Dot>> &Dotss, int n)
    {
        for (int s = 0; s < Dotss.size(); ++s)
            dot_initial(Dotss[s], n, s);
    }

    void dot_initial(std::deque<Dot> &Dots, int n, int s)
    {
        if (Dots.size() < O.size() + 1)
            Dots.resize(O.size() + 1);
        for (int w = 0; w <= O.size(); ++w)
        {
            Dots[w].n = n;
            Dots[w].s = s;
            Dots[w].w = w;
            Dots[w].s_sz = 0;
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
        Q.n = -3;
        Q.s = 0;
        Q.w = 0;
    }

    void renew_source(double newval, Dot *newadr, double &maxval, std::deque<Dot *> &maxadr)
    {
        if (newval >= maxval)
        {
            if (newval > maxval)
            {
                maxval = newval;
                maxadr.clear();
            }
            maxadr.emplace_back(newadr);
        }
    }

    void renew_source(double newval, int64_t sr1, int32_t lambda, double &maxval, std::deque<std::pair<int64_t, int32_t>> &maxadr)
    {
        if (newval >= maxval)
        {
            if (newval > maxval)
            {
                maxval = newval;
                maxadr.clear();
            }
            maxadr.emplace_back(sr1, lambda);
        }
    }

    void Mix()
    {
        if (O.size() + 1 > crossglobaldata.Ethres.size())
        {
            crossglobaldata.Ethres.resize(O.size() + 1);
            crossglobaldata.C.resize(O.size() + 1);
            for (auto &edge : global_circuits)
                edge.Cdelta.resize(O.size() + 1);
            for (auto &edge : global_crosses)
                edge.Cdelta.resize(O.size() + 1);
            for (int n = 0; n < nodes.size(); ++n)
            {
                nodes[n].AdeltaDot.resize(O.size() + 1);
                nodes[n].AdeltaCross.resize(O.size() + 1);
                nodes[n].AdeltaCircuit.resize(O.size() + 1);
            }
        }
        for (int n = 0; n < nodes.size(); ++n)
        {
            dot_initial(nodes[n].A, -n - 4);
            dot_initial(nodes[n].B, -1, 0);
            for (int w = 0; w <= O.size(); ++w)
                nodes[n].A[0][w].val = -inf;
        }
        for (auto &scc : sccs)
        {
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
                alloc_initial(edge->D0, edge->n, scc.nodes.size() - 1);
            for (size_t w = 0; w <= O.size(); ++w)
            {
                for (auto node : scc.nodes)
                {
                    if (w == 0)
                        source_max(node->B[w], {}, {});
                    else
                        source_max(node->B[w], {&node->B[w - 1], &node->A[scc.nodes.size() - 1][w - 1]}, {node->B[w - 1].val + node->ue, node->A[scc.nodes.size() - 1][w - 1].val + node->ve});

                    double oldval = node->A[0][w].val;
                    renew_source(node->B[w].val, &node->B[w], node->A[0][w].val, node->AdeltaDot[w]);
                    if (w == 0)
                        renew_source(node->Abar.val, &node->Abar, node->A[0][w].val, node->AdeltaDot[w]);
                    if (node->A[0][w].val > oldval)
                    {
                        for (auto cross : node->AdeltaCross[w])
                            cross->Cdelta[w].clear();
                        node->AdeltaCross[w].clear();
                    }
                }
                for (auto edge : scc.local_circuits)
                    CircuitIteration(w, edge, scc.nodes.size());
                for (auto edge : scc.global_circuits)
                    CircuitIterationGlobal(w, edge, scc.nodes.size());
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
                CrossIterationGlobal(edge, scc.nodes.size());
        }
        std::vector<Dot *> adrs;
        std::vector<double> vals;
        for (size_t n = 0; n < targets.size(); ++n)
        {
            adrs.push_back(&targets[n]->A.back()[O.size()]);
            vals.push_back(targets[n]->A.back()[O.size()].val);
        }
        source_max(Q, adrs.begin(), vals.begin(), vals.end());
    }

    void CrossIteration(EdgeLocalCross *edge, int scc_sz)
    {
        std::deque<Dot> &A = edge->tail->A.back();
        Dot **E = edge->E;
        Dot **F = edge->F;
        Dot **G = edge->G;
        std::deque<Dot> &hA = edge->head->A.front();
        std::deque<std::deque<Dot *>> &AdeltaDot = edge->head->AdeltaDot;
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
                    source_max(G[s][w], {&F[s][w], &E[s][w], &A[w]}, {F[s][w].val, E[s][w].val, A[w].val + (s == 0 ? 0 : edge->vfp + (s - 1) * edge->ufp) + edge->T});
                else
                    source_max(G[s][w], {&F[s][w], &E[s][w], &G[s - 1][w - 1], &A[w]}, {F[s][w].val, E[s][w].val, G[s - 1][w - 1].val + edge->gamma[BroWheel::base2int(edge->nameseq.seq[s - 1])][BroWheel::base2int(O[w - 1])], A[w].val + (s == 0 ? 0 : edge->vfp + (s - 1) * edge->ufp) + edge->T});
                renew_source(G[s][w].val, &G[s][w], hA[w].val, AdeltaDot[w]);
            }
    }

    void CrossIterationGlobal(EdgeGlobalCross *edge, int scc_sz)
    {
        std::deque<Dot> &A = edge->tail->A.back();
        auto &simnodes = crossglobaldata.simnodes;
        auto &E0 = crossglobaldata.E0;
        auto &E = crossglobaldata.E;
        std::deque<double> &C = crossglobaldata.C, &Ethres = crossglobaldata.Ethres;
        std::deque<std::deque<std::pair<int64_t, int>>> &Cdelta = edge->Cdelta;

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
            C[w] = edge->head->A[0][w].val;
            renew_source(simnodes[lambda].FG[w].G, simnodes[lambda].sr1, lambda, C[w], Cdelta[w]);
        }

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
                    renew_source(G_val, simnodes[lambda].sr1, lambda, C[w], Cdelta[w]);
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

        for (int w = 0; w <= O.size(); ++w)
            if (!Cdelta[w].empty())
            {
                if (C[w] > edge->head->A[0][w].val)
                {
                    edge->head->A[0][w].val = C[w];
                    edge->head->AdeltaDot[w].clear();
                    for (auto cross : edge->head->AdeltaCross[w])
                        cross->Cdelta[w].clear();
                    edge->head->AdeltaCross[w].clear();
                }
                edge->head->AdeltaCross[w].emplace_back(edge);
            }
    }

    void CircuitIteration(size_t w, EdgeLocalCircuit *edge, int scc_sz)
    {
        std::deque<Dot> &A = edge->tail->A.back();
        Dot *D = edge->D;
        Dot **E = edge->E;
        Dot **F0 = edge->F0;
        Dot **G0 = edge->G0;
        Dot **G = edge->G;
        double C = edge->head->A[0][w].val;
        std::deque<Dot *> &AdeltaDot = edge->head->AdeltaDot[w];

        if (w > 0)
            source_max(D[w - 1], {&A[w - 1]}, {A[w - 1].val + edge->T});
        for (size_t s = 0; s <= edge->nameseq.seq.size(); ++s)
        {
            if (w > 0)
            {
                if (s > 0)
                    source_max(G[s][w - 1], {&G0[s][w - 1], &D[w - 1], &A[w - 1]}, {G0[s][w - 1].val, D[w - 1].val + edge->vf + (s - 1) * edge->uf, A[w - 1].val + edge->vfp + (s - 1) * edge->ufp + edge->T});
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

            renew_source(G0[s][w].val, &G0[s][w], C, AdeltaDot);
        }

        if (C > edge->head->A[0][w].val)
        {
            edge->head->A[0][w].val = C;
            for (auto cross : edge->head->AdeltaCross[w])
                cross->Cdelta[w].clear();
            edge->head->AdeltaCross[w].clear();
        }
    }

    void CircuitIterationGlobal(size_t w, EdgeGlobalCircuit *edge, int scc_sz)
    {
        if (w == 0)
            return;

        std::deque<Dot> &A = edge->tail->A.back(), &B = edge->tail->B;
        std::deque<SNC> &sncs = edge->sncs;
        std::deque<std::pair<int64_t, int>> &Cdelta = edge->Cdelta[w];

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
        double C = edge->head->A[0][w].val;
        renew_source(sncs.front().G0, sncs.front().sr1, sncs.front().lambda, C, Cdelta);

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
            else
                renew_source(jump->G0, jump->sr1, jump->lambda, C, Cdelta);
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

        if (!Cdelta.empty())
        {
            if (C > edge->head->A[0][w].val)
            {
                edge->head->A[0][w].val = C;
                edge->head->AdeltaDot[w].clear();
                for (auto cross : edge->head->AdeltaCross[w])
                    cross->Cdelta[w].clear();
                edge->head->AdeltaCross[w].clear();
                for (auto circuit : edge->head->AdeltaCircuit[w])
                    circuit->Cdelta[w].clear();
                edge->head->AdeltaCircuit[w].clear();
            }
            edge->head->AdeltaCircuit[w].emplace_back(edge);
        }
    }

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
                for (auto edge : node.AdeltaCross[w])
                    edge->Cdelta[w].clear();
                node.AdeltaCross[w].clear();
                for (auto edge : node.AdeltaCircuit[w])
                    edge->Cdelta[w].clear();
                node.AdeltaCircuit[w].clear();
                node.AdeltaDot[w].clear();
            }
        clear();
    }

    void GlobalTrack(Dot *M)
    {
        for (auto &edge : nodes[-M->n - 4].AdeltaCross[M->w])
        {
            std::deque<Dot> &A = edge->tail->A.back();
            auto &tracknodes = edge->tracknodes;
            for (auto &suf : edge->Cdelta[M->w])
            {
                size_t start = edge->browheel.SimSuffix(suf.first);
                if (tracknodes.empty())
                    tracknodes.emplace_back(this, (TrackNode *)NULL, (TrackNode *)NULL, edge->n, edge->browheel.sequence.size() - 1 - start - suf.second, O.size(), 0);
                auto tracknode = &tracknodes.front();
                for (int w = tracknode->tau; w <= M->w; ++w)
                {
                    if (w == 0)
                        source_max(tracknode->E[w], {}, {});
                    else
                        source_max(tracknode->E[w], {&tracknode->E[w - 1], &tracknode->G[w - 1]}, {tracknode->E[w - 1].val + edge->ue, tracknode->G[w - 1].val + edge->ve});
                    source_max(tracknode->F[w], {}, {});
                    source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &A[w]}, {tracknode->E[w].val, tracknode->F[w].val, A[w].val + edge->T});
                }
                tracknode->tau = std::max(tracknode->tau, M->w + 1);
                for (int i = suf.second - 1; i >= 0; --i)
                {
                    int c = edge->browheel.sequence(start + i);
                    if (!tracknode->itcs[c - 2])
                    {
                        tracknodes.emplace_back(this, tracknode, (TrackNode *)NULL, edge->n, edge->browheel.sequence.size() - 1 - start - i, O.size(), tracknode->lambda + 1);
                        tracknode->itcs[c - 2] = &tracknodes.back();
                    }
                    tracknode = tracknode->itcs[c - 2];
                    for (int w = tracknode->tau; w <= M->w; ++w)
                    {
                        if (w == 0)
                            source_max(tracknode->E[w], {}, {});
                        else
                            source_max(tracknode->E[w], {&tracknode->E[w - 1], &tracknode->G[w - 1]}, {tracknode->E[w - 1].val + edge->ue, tracknode->G[w - 1].val + edge->ve});
                        source_max(tracknode->F[w], {&tracknode->itp->F[w], &tracknode->itp->G[w]}, {tracknode->itp->F[w].val + edge->uf, tracknode->itp->G[w].val + edge->vf});
                        if (w == 0)
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w]}, {tracknode->E[w].val, tracknode->F[w].val});
                        else
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &tracknode->itp->G[w - 1]}, {tracknode->E[w].val, tracknode->F[w].val, tracknode->itp->G[w - 1].val + edge->gamma[c][BroWheel::base2int(O[w - 1])]});
                    }
                    tracknode->tau = std::max(tracknode->tau, M->w + 1);
                }
                nodes[-M->n - 4].AdeltaDot[M->w].emplace_back(&tracknode->G[M->w]);
            }
        }
        for (auto &edge : nodes[-M->n - 4].AdeltaCircuit[M->w])
        {
            std::deque<Dot> &A = edge->tail->A.back();
            auto &tracknodes = edge->tracknodes;
            for (auto &suf : edge->Cdelta[M->w])
            {
                size_t start = edge->browheel.SimSuffix(suf.first);
                if (tracknodes.empty())
                {
                    alloc_initial(edge->G, edge->n, edge->browheel.sequence.size() - 1 - start - suf.second);
                    for (size_t w = 0; w <= O.size(); ++w)
                        edge->G[w].lambda = 0;
                    tracknodes.emplace_back(this, (TrackNode *)NULL, (TrackNode *)NULL, edge->n, edge->browheel.sequence.size() - 1 - start - suf.second, O.size(), 0);
                }
                auto tracknode = &tracknodes.front();
                for (int w = tracknode->tau; w <= M->w; ++w)
                {
                    if (w > 0)
                        source_max(edge->G[w - 1], {&tracknode->G[w - 1], &A[w - 1]}, {tracknode->G[w - 1].val, A[w - 1].val + edge->T});
                    if (w == 0)
                        source_max(tracknode->E[w], {}, {});
                    else
                        source_max(tracknode->E[w], {&tracknode->E[w - 1], &edge->G[w - 1]}, {tracknode->E[w - 1].val + edge->ue, edge->G[w - 1].val + edge->ve});
                    source_max(tracknode->F[w], {}, {});
                    source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w]}, {tracknode->E[w].val, tracknode->F[w].val});
                }
                tracknode->tau = std::max(tracknode->tau, M->w + 1);
                for (int i = suf.second - 1; i >= 0; --i)
                {
                    int c = edge->browheel.sequence(start + i);
                    if (!tracknode->itcs[c - 2])
                    {
                        tracknodes.emplace_back(this, tracknode, (TrackNode *)NULL, edge->n, edge->browheel.sequence.size() - 1 - start - i, O.size(), tracknode->lambda + 1);
                        tracknode->itcs[c - 2] = &tracknodes.back();
                    }
                    tracknode = tracknode->itcs[c - 2];
                    for (int w = tracknode->tau; w <= M->w; ++w)
                    {
                        if (w == 0)
                            source_max(tracknode->E[w], {}, {});
                        else
                            source_max(tracknode->E[w], {&tracknode->E[w - 1], &tracknode->G[w - 1]}, {tracknode->E[w - 1].val + edge->ue, tracknode->G[w - 1].val + edge->ve});
                        source_max(tracknode->F[w], {&tracknode->itp->F[w], &tracknode->itp->G[w]}, {tracknode->itp->F[w].val + edge->uf, tracknode->itp->G[w].val + edge->vf});
                        if (w == 0)
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w]}, {tracknode->E[w].val, tracknode->F[w].val});
                        else if (tracknode->itp == &tracknodes.front())
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &edge->G[w - 1]}, {tracknode->E[w].val, tracknode->F[w].val, edge->G[w - 1].val + edge->gamma[c][BroWheel::base2int(O[w - 1])]});
                        else
                            source_max(tracknode->G[w], {&tracknode->E[w], &tracknode->F[w], &tracknode->itp->G[w - 1]}, {tracknode->E[w].val, tracknode->F[w].val, tracknode->itp->G[w - 1].val + edge->gamma[c][BroWheel::base2int(O[w - 1])]});
                    }
                    tracknode->tau = std::max(tracknode->tau, M->w + 1);
                }
                nodes[-M->n - 4].AdeltaDot[M->w].emplace_back(&tracknode->G[M->w]);
            }
        }
        M->s_sz = nodes[-M->n - 4].AdeltaDot[M->w].size();
        M->sources = heap_alloc<Dot *>(M->s_sz);
        for (int i = 0; i < M->s_sz; ++i)
            M->sources[i] = nodes[-M->n - 4].AdeltaDot[M->w][i];
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
