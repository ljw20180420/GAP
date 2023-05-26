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

struct Doo
{
    int w;
    double F;
    double G;
};

struct SimNode
{
    Memory *memory;
    int c;
    size_t sr1;
    size_t sr2;
    Doo *FG;
    int sz;
    char *now;
    size_t remain;
    size_t chunk_id;

    SimNode(Memory *memory_, int c_, size_t sr1_, size_t sr2_, int sz_, int alloc)
        : memory(memory_)
    {
        c = c_;
        sr1 = sr1_;
        sr2 = sr2_;

        now = memory->now;
        remain = memory->remain;
        chunk_id = memory->chunk_id;
        sz = sz_;
        FG = memory->heap_alloc<Doo>(alloc);
    }

    ~SimNode()
    {
        memory->now = now;
        memory->remain = remain;
        memory->chunk_id = chunk_id;
    }
};

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
    {
        E = E_;
        F0 = F0_;
        G0 = G0_;
        hatG = hatG_;
        c = c_;
        idl = idl_;
        lambda = lambda_;
        sr1 = sr1_;
        sr2 = sr2_;
        cid = cid_;
        itp = itp_;
        jump = jump_;
    }
};

struct Align : Memory, Graph
{
    const static int ww = 3;

    std::string O;
    std::string Oname;
    Dot Q;
    std::vector<int> target_scc_szes;
    std::vector<int> node_scc_szes;
    std::ofstream fout;
    size_t max_id=0;

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
        clear();
        for (auto &scc : sccs)
        {
            for (auto node : scc.nodes)
            {
                alloc_initial(node->A, nodes.data() - node - 4, scc.nodes.size() - 1);
                alloc_initial(node->B, -1, 0);
                node->AdeltaCross.resize(O.size() + 1);
                node->AdeltaCircuit.resize(O.size() + 1);
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
                edge->Cdelta.resize(O.size() + 1);
            }
            for (size_t w = 0; w <= O.size(); ++w)
            {
                for (auto edge : scc.local_circuits)
                    CircuitIteration(w, edge, scc.nodes.size());
                for (auto edge : scc.global_circuits)
                    CircuitIterationGlobal(w, edge, scc.nodes.size());
                for (auto node : scc.nodes)
                {
                    if (w == 0)
                        source_max(node->B[w], {}, {});
                    else
                        source_max(node->B[w], {&node->B[w - 1], &node->A[scc.nodes.size() - 1][w - 1]}, {node->B[w - 1].val + node->ue, node->A[scc.nodes.size() - 1][w - 1].val + node->ve});
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
                    node->AdeltaCross[w].clear();
                    for (auto edge : node->in_global_crosses)
                        if (edge->C[w] >= node->A[0][w].val)
                        {
                            if (edge->C[w] > node->A[0][w].val)
                            {
                                node->A[0][w].val = edge->C[w];
                                node->A[0][w].s_sz = 0;
                                node->AdeltaCross[w].clear();
                            }
                            node->AdeltaCross[w].insert_or_assign(edge, std::move(edge->Cdelta[w]));
                        }
                    node->AdeltaCircuit[w].clear();
                    for (auto edge : node->in_global_circuits)
                        if (edge->C[w] >= node->A[0][w].val)
                        {
                            if (edge->C[w] > node->A[0][w].val)
                            {
                                node->A[0][w].val = edge->C[w];
                                node->A[0][w].s_sz = 0;
                                node->AdeltaCircuit[w].clear();
                            }
                            node->AdeltaCircuit[w].insert_or_assign(edge, std::move(edge->Cdelta[w]));
                        }
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
                edge->Cdelta.resize(O.size() + 1);
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
        Dot **A = edge->tail->A;
        std::deque<SimNode> simnodes;
        simnodes.emplace_front(this, -1, 0, edge->browheel.sequence.size() - 1, O.size() + 1, O.size() + 1);
        auto &simzero = simnodes.back();
        std::vector<Do> E0;
        for (size_t w = 0; w <= O.size(); ++w)
        {
            E0.emplace_back(w, w > 0 ? std::max(E0[w - 1].E + edge->ue, simzero.FG[w - 1].G + edge->ve) : -inf);
            simzero.FG[w].w = w;
            simzero.FG[w].F = -inf;
            simzero.FG[w].G = std::max(std::max(E0[w].E, simzero.FG[w].F), A[scc_sz - 1][w].val + edge->T);
            edge->C[w] = simzero.FG[w].G;
            edge->Cdelta[w].clear();
            edge->Cdelta[w].emplace_back(simzero.sr1, simnodes.size() - 1);
        }
        std::deque<SimNode>::iterator simprev = simnodes.begin(), simnode;
        int c = 1;
        int64_t sr1, sr2;
        do
        {
            ++c;
            sr1 = simprev->sr1;
            sr2 = simprev->sr2;
            edge->browheel.PreRange(sr1, sr2, c);
        } while (sr1 > sr2 && c + 1 < 7);
        if (sr1 > sr2)
            simnodes.clear();
        else
        {
            simnodes.emplace_front(this, c, sr1, sr2, simprev->sz, O.size() + 1 - simprev->FG[0].w);
            simnode = simnodes.begin();
        }
        while (!simnodes.empty())
        {
            std::vector<Do> E;
            for (size_t i = 0, j = 0, w = simprev->FG[j].w;;)
            {
                double E_val = std::max((!E.empty() && size_t(E.back().w) == w - 1) ? E.back().E + edge->ue : -inf, (i > 0 && size_t(simnode->FG[i - 1].w) == w - 1) ? simnode->FG[i - 1].G + edge->ve : -inf);
                double G_val = j > 0 && size_t(simprev->FG[j - 1].w) == w - 1 ? simprev->FG[j - 1].G + edge->gamma[simnode->c][BroWheel::base2int(O[w - 1])] : -inf;
                if (E_val >= std::max(E0[w].E, simzero.FG[w].G + (O.size() - w) * edge->tail->ve))
                {
                    E.emplace_back(w, E_val);
                    G_val = E_val > G_val ? E_val : G_val;
                }
                double F_val = j < size_t(simprev->sz) && size_t(simprev->FG[j].w) == w ? std::max(simprev->FG[j].F + edge->uf, simprev->FG[j].G + edge->vf) : -inf;
                G_val = F_val > G_val ? F_val : G_val;
                if (G_val >= simzero.FG[w].G)
                {
                    simnode->FG[i].w = w;
                    simnode->FG[i].F = F_val;
                    simnode->FG[i].G = G_val;
                    if (G_val >= edge->C[w])
                    {
                        if (G_val > edge->C[w])
                        {
                            edge->C[w] = G_val;
                            edge->Cdelta[w].clear();
                        }
                        edge->Cdelta[w].emplace_back(simnode->sr1, simnodes.size() - 1);
                    }
                    ++i;
                }
                while (size_t(simprev->FG[j].w) <= w && j + 1 < size_t(simprev->sz))
                    ++j;
                if (w < O.size() && (size_t(simprev->FG[j - 1].w) == w || size_t(simprev->FG[j].w) == w || (i > 0 && size_t(simnode->FG[i - 1].w) == w) || (!E.empty() && size_t(E.back().w) == w)))
                    ++w;
                else if (size_t(simprev->FG[j].w) > w)
                    w = simprev->FG[j].w;
                else
                {
                    simnode->sz = i;
                    break;
                }
            }
            do
            {
                if (simnode->sr1 <= simnode->sr2 && simnode->sz > 0)
                    c = 1;
                else
                {
                    do
                    {
                        c = simnodes.front().c;
                        simnodes.pop_front();
                    } while (c == 6 && !simnodes.empty());
                    if (simnodes.empty())
                        break;
                }
                simprev = simnodes.begin();
                do
                {
                    ++c;
                    sr1 = simprev->sr1;
                    sr2 = simprev->sr2;
                    edge->browheel.PreRange(sr1, sr2, c);
                } while (sr1 > sr2 && c + 1 < 7);
                if (sr1 <= sr2)
                    simnodes.emplace_front(this, c, sr1, sr2, simprev->sz, O.size() + 1 - simprev->FG[0].w);
                else
                    simprev->sz=0;
                simnode = simnodes.begin();
            } while (sr1 > sr2);
        }
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
        Dot **A = edge->tail->A;
        if (w == 0)
        {
            edge->sncs.emplace_back(-inf, -inf, -inf, -inf, 0, 5, 0, 0, edge->browheel.sequence.size() - 1, 1, (SNC *)NULL, (SNC *)NULL);
            for (int c = 2; c <= 6; ++c)
            {
                int64_t sr1 = edge->sncs.front().sr1;
                int64_t sr2 = edge->sncs.front().sr2;
                edge->browheel.PreRange(sr1, sr2, c);
                edge->sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, 1, sr1, sr2, 0, &edge->sncs.front(), edge->sncs.front().jump);
                edge->sncs.front().jump=&edge->sncs.back();
            }
            edge->C[0] = -inf;
            edge->Cdelta[0].clear();
            return;
        }
        edge->sncs.front().hatG = std::max(edge->sncs.front().G0, A[scc_sz - 1][w - 1].val + edge->T);
        for (auto prejump = &edge->sncs.front(), jump = prejump->jump; jump; jump = prejump->jump)
        {
            jump->hatG = jump->G0;
            if (jump->itp != &edge->sncs.front() && jump->itp->hatG == -inf && jump->hatG == -inf)
                prejump->jump=jump->jump;
            else
                prejump=jump;
        }
        edge->sncs.front().E = std::max(edge->sncs.front().hatG + edge->ve, edge->sncs.front().E + edge->ue);
        edge->sncs.front().F0 = -inf;
        edge->sncs.front().G0 = std::max(edge->sncs.front().E, edge->sncs.front().F0);
        edge->C[w] = edge->sncs.front().G0;
        edge->Cdelta[w].clear();
        edge->Cdelta[w].emplace_back(edge->sncs.front().sr1, 0);

        for (auto jump = edge->sncs.front().jump; jump; jump = jump->jump)
        {
            jump->E = std::max(jump->hatG + edge->ve, jump->E + edge->ue);
            jump->F0 = std::max(jump->itp->G0 + edge->vf, jump->itp->F0 + edge->uf);
            jump->G0 = std::max(std::max(jump->E, jump->F0), jump->itp->hatG + edge->gamma[jump->c][BroWheel::base2int(O[w - 1])]);
            if (jump->G0 < edge->sncs.front().G0)
            {
                jump->E = -inf;
                jump->F0 = -inf;
                jump->G0 = -inf;
            }
            if (jump->G0 >= edge->C[w])
            {
                if (jump->G0 > edge->C[w])
                {
                    edge->C[w] = jump->G0;
                    edge->Cdelta[w].clear();
                }
                edge->Cdelta[w].emplace_back(jump->sr1, jump->lambda);
            }
            if (jump->G0 != -inf && jump->hatG == -inf)
            {
                if (jump->cid == 0)
                {
                    jump->cid=edge->sncs.size();
                    for (int c = 2; c <= 6; ++c)
                    {
                        int64_t sr1 = jump->sr1;
                        int64_t sr2 = jump->sr2;
                        edge->browheel.PreRange(sr1, sr2, c);
                        if (sr1 <= sr2)
                        {
                            edge->sncs.emplace_back(-inf, -inf, -inf, -inf, c, 0, jump->lambda+1, sr1, sr2, 0, jump, (SNC *)NULL);
                            ++jump->idl;
                        }
                    }   
                }
                for (int idi=0; idi<jump->idl; ++idi)
                    if (edge->sncs[jump->cid+idi].hatG == -inf)
                    {
                        edge->sncs[jump->cid+idi].jump=jump->jump;
                        jump->jump=&edge->sncs[jump->cid+idi];
                    }
            }
        }
        if (w == O.size())
            edge->sncs.clear();
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
                    tracknodes.emplace_back(this, (TrackNode *)NULL, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size()-1-start-suf.second, O.size(), 0);
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
                        tracknodes.emplace_back(this, tracknode, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size()-1-start-i, O.size(), tracknode->lambda + 1);
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
                    alloc_initial(pair.first->G, pair.first->n, pair.first->browheel.sequence.size()-1-start-suf.second);
                    for (size_t w=0; w<=O.size(); ++w)
                        pair.first->G[w].lambda=0;
                    tracknodes.emplace_back(this, (TrackNode *)NULL, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size()-1-start-suf.second, O.size(), 0);
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
                        tracknodes.emplace_back(this, tracknode, (TrackNode *)NULL, pair.first->n, pair.first->browheel.sequence.size()-1-start-i, O.size(), tracknode->lambda + 1);
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
