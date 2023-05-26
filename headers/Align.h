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
            int shiftFG;

            SimNode(int c_, int64_t sr1_, int64_t sr2_, int shiftFG_)
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

            Do(int w_, double E_)
            : w(w_), E(E_)
            {
            }
        };

        std::deque<Do> E0, E;
        std::deque<double> Ethres;
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
            ptr[w].visit = false;
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
            Dots[w].visit = false;
        }
    }

    Align(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel, size_t chunk_sz_, std::string file) : Graph(argc, argv, file2seq, file2browheel), fout(file, std::ifstream::binary)
    {
        Initial(chunk_sz_);
        for (auto &node : nodes)
        {
            node.Abar.visit = false;
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
            for (int n = 0; n < nodes.size(); ++n)
            {
                nodes[n].AdeltaDot.resize(O.size() + 1);
                nodes[n].AdeltaGlobal.resize(O.size() + 1);
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

                    node->updateA0(w, node->B[w].val, &node->B[w]);
                    if (w == 0)
                        node->updateA0(w, node->Abar.val, &node->Abar);
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
                edge->head->updateA0(w, G[s][w].val, &G[s][w]);
            }
    }

    void CrossIterationGlobal(EdgeGlobalCross *edge, int scc_sz)
    {
        std::deque<Dot> &A = edge->tail->A.back();
        std::deque<CrossGlobalData::SimNode> &simnodes = crossglobaldata.simnodes;
        std::deque<CrossGlobalData::Do> &E0 = crossglobaldata.E0, &E = crossglobaldata.E;
        std::deque<double> &Ethres = crossglobaldata.Ethres;
        std::deque<CrossGlobalData::Doo> &FG = crossglobaldata.FG;

        int c = -1;
        int64_t sr1 = 0, sr2 = edge->browheel.sequence.size() - 1;
        crossglobaldata.emplace_back(c, sr1, sr2);
        for (size_t w = 0; w <= O.size(); ++w)
        {
            double tgw = (O.size() > w ? (O.size() - w - 1) * edge->tail->ue + edge->tail->ve : 0);
            E0.emplace_back(w, w > 0 ? std::max(E0[w-1].E + edge->ue, FG[w - 1].G + edge->ve) : -inf);
            if (E0[w].E + edge->ue < std::max(E0[w].E, A[w].val + edge->T) + edge->tail->ve && E0[w].E + (O.size() - w) * edge->ue < std::max(E0[w].E, A[w].val + edge->T) + tgw)
                E0[w].E = -inf;
            FG.emplace_back(w, -inf, std::max(E0[w].E, A[w].val + edge->T));
            Ethres[w] = std::max(E0[w].E, FG[w].G + std::min(0.0, std::min(edge->tail->ve - edge->ue, tgw - (O.size() - w) * edge->ue)));
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
                    edge->browheel.PreRange(sr1, sr2, c);
                } while (sr1 > sr2 && c < 6);
                if (sr1 <= sr2)
                {
                    crossglobaldata.emplace_back(c, sr1, sr2);
                    break;
                }
                else if(FG.size() > simnodes.back().shiftFG)
                    FG.erase(FG.begin() + simnodes.back().shiftFG, FG.end());
            } while (true);
            if (simnodes.empty())
                break;

            int64_t idx = simnodes[simnodes.size() - 2].shiftFG;
            for (int j = 0, w = FG[idx].w;;)
            {
                E.emplace_back(w, std::max((!E.empty() && E.back().w == w - 1) ? E.back().E + edge->ue : -inf, (simnodes.back().shiftFG < FG.size() && FG.back().w == w - 1) ? FG.back().G + edge->ve : -inf));
                if (E.back().E < Ethres[w])
                    E.back().E = -inf;

                double F_val = idx < simnodes.back().shiftFG && FG[idx].w == w ? std::max(FG[idx].F + edge->uf, FG[idx].G + edge->vf) : -inf;
                double G_val = j > 0 && FG[idx - 1].w == w - 1 ? FG[idx - 1].G + edge->gamma[simnodes.back().c][BroWheel::base2int(O[w - 1])] : -inf;
                G_val = std::max(std::max(G_val, E.back().E), F_val);
                if (G_val >= FG[w].G)
                {
                    FG.emplace_back(w, F_val, G_val);
                    edge->head->updateA0(w, G_val, edge, simnodes.back().sr1, simnodes.size() - 1);
                }
                if (idx < simnodes.back().shiftFG && FG[idx].w <= w)
                {
                    ++j;
                    ++idx;
                }
                if (w < O.size() && (FG[idx - 1].w == w || (simnodes.back().shiftFG < FG.size() && FG.back().w == w) || (!E.empty() && E.back().w == w && E.back().E > -inf)))
                    ++w;
                else if (idx < simnodes.back().shiftFG)
                    w = FG[idx].w;
                else
                    break;
            }

            E.clear();
        }

        E0.clear();
    }

    void CircuitIteration(size_t w, EdgeLocalCircuit *edge, int scc_sz)
    {
        std::deque<Dot> &A = edge->tail->A.back();
        Dot *D = edge->D;
        Dot **E = edge->E;
        Dot **F0 = edge->F0;
        Dot **G0 = edge->G0;
        Dot **G = edge->G;

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
            edge->head->updateA0(w, G0[s][w].val, &G0[s][w]);
        }
    }

    void CircuitIterationGlobal(size_t w, EdgeGlobalCircuit *edge, int scc_sz)
    {
        if (w == 0)
            return;

        std::deque<Dot> &A = edge->tail->A.back(), &B = edge->tail->B;
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
        edge->head->updateA0(w, sncs.front().G0, edge, sncs.front().sr1, sncs.front().lambda);

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
                edge->head->updateA0(w, jump->G0, edge, jump->sr1, jump->lambda);
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
    }

    void GetMinimalGraph()
    {
        int Onsz = Oname.size();
        fout.write((char *)&Onsz, sizeof(Onsz));
        fout.write(Oname.data(), sizeof(char) * Onsz);
        std::queue<Dot *> dots;
        std::deque<Dot *> visits;
        size_t id = 0;
        Q.id = id;
        Q.visit = true;
        dots.push(&Q);
        while (!dots.empty())
        {
            Dot *M = dots.front();
            dots.pop();
            visits.push_back(M);
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
                if (!M->sources[i]->visit)
                {
                    M->sources[i]->id = ++id;
                    M->sources[i]->visit = true;
                    dots.push(M->sources[i]);
                }
                fout.write((char *)&M->sources[i]->id, sizeof(M->sources[i]->id));
            }
        }
        if (id > max_id)
            max_id = id;

        while(!visits.empty())
        {
            visits.back()->visit = false;
            visits.pop_back();
        }
        for (auto &edge : global_crosses)
            edge.tracktree.clear();
        for (auto &edge : global_circuits)
            edge.tracktree.clear();
        for (auto &node : nodes)
            node.clearAdelta(O.size());
        clear();
    }

    void GlobalTrack(Dot *M)
    {
        for (auto &globalsuffix : nodes[-M->n - 4].AdeltaGlobal[M->w])
        {
            EdgeGlobal *edge = globalsuffix.edge;
            std::deque<Dot> &A = edge->tail->A.back();
            TrackTree &tracktree = edge->tracktree;
            std::deque<TrackTree::TrackNode> &tracknodes = tracktree.tracknodes;
            std::deque<Dot> &dots = tracktree.dots;
            size_t start = edge->browheel.SimSuffix(globalsuffix.start);
            if (tracktree.tracknodes.empty())
            {
                tracktree.W = O.size();
                tracktree.emplace_back(-1, edge->n, edge->browheel.sequence.size() - 1 - start - globalsuffix.lambda, 0);
            }
            tracktree.setidx(0);
            if (edge->head->scc_id != edge->tail->scc_id)
            {
                for (int w = tracknodes[tracktree.idx].tau; w <= M->w; ++w)
                {
                    if (w == 0)
                        source_max(dots[tracktree.shiftE + w], {}, {});
                    else
                        source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                    source_max(dots[tracktree.shiftF + w], {}, {});
                    source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], &A[w]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, A[w].val + edge->T});
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, M->w + 1);
                for (int i = globalsuffix.lambda - 1; i >= 0; --i)
                {
                    int c = edge->browheel.sequence(start + i);
                    if (tracknodes[tracktree.idx].cidxs[c - 2] < 0)
                    {
                        tracktree.emplace_back(-1, edge->n, edge->browheel.sequence.size() - 1 - start - i, dots[tracktree.shiftE].lambda + 1);
                        tracknodes[tracktree.idx].cidxs[c - 2] = tracknodes.size() - 1;
                    }
                    tracktree.setidx(tracknodes[tracktree.idx].cidxs[c - 2]);
                    for (int w = tracknodes[tracktree.idx].tau; w <= M->w; ++w)
                    {
                        if (w == 0)
                            source_max(dots[tracktree.shiftE + w], {}, {});
                        else
                            source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                        source_max(dots[tracktree.shiftF + w], {&dots[tracktree.shiftPF + w], &dots[tracktree.shiftPG + w]}, {dots[tracktree.shiftPF + w].val + edge->uf, dots[tracktree.shiftPG + w].val + edge->vf});
                        if (w == 0)
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val});
                        else
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], &dots[tracktree.shiftPG + w - 1]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, dots[tracktree.shiftPG + w - 1].val + edge->gamma[c][BroWheel::base2int(O[w - 1])]});
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, M->w + 1);
                }
                nodes[-M->n - 4].AdeltaDot[M->w].emplace_back(&dots[tracktree.shiftG + M->w]);
            }
            else
            {
                for (int w = tracknodes[tracktree.idx].tau; w <= M->w; ++w)
                {
                    if (w > 0)
                        source_max(dots[tracktree.shiftG + w - 1], {&dots[tracktree.shiftE + w - 1], &A[w - 1]}, {dots[tracktree.shiftE + w - 1].val, A[w - 1].val + edge->T});
                    if (w == 0)
                        source_max(dots[tracktree.shiftE + w], {}, {});
                    else
                        source_max(dots[tracktree.shiftE + w], {&dots[tracktree.shiftE + w - 1], &dots[tracktree.shiftG + w - 1]}, {dots[tracktree.shiftE + w - 1].val + edge->ue, dots[tracktree.shiftG + w - 1].val + edge->ve});
                    source_max(dots[tracktree.shiftF + w], {}, {});
                }
                tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, M->w + 1);
                for (int i = globalsuffix.lambda - 1; i >= 0; --i)
                {
                    int c = edge->browheel.sequence(start + i);
                    if (tracknodes[tracktree.idx].cidxs[c - 2] < 0)
                    {
                        tracktree.emplace_back(-1, edge->n, edge->browheel.sequence.size() - 1 - start - i, dots[tracktree.shiftE].lambda + 1);
                        tracknodes[tracktree.idx].cidxs[c - 2] = tracknodes.size() - 1;
                    }
                    tracktree.setidx(tracknodes[tracktree.idx].cidxs[c - 2]);
                    for (int w = tracknodes[tracktree.idx].tau; w <= M->w; ++w)
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
                            source_max(dots[tracktree.shiftG + w], {&dots[tracktree.shiftE + w], &dots[tracktree.shiftF + w], &dots[tracktree.shiftPG + w - 1]}, {dots[tracktree.shiftE + w].val, dots[tracktree.shiftF + w].val, dots[tracktree.shiftPG + w - 1].val + edge->gamma[c][BroWheel::base2int(O[w - 1])]});
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, M->w + 1);
                }
                nodes[-M->n - 4].AdeltaDot[M->w].emplace_back(&dots[tracktree.shiftG + M->w]);
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
