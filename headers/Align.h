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
#include "Memory.h"

struct Do
{
    int w;
    double val;
};

template <typename U>
struct SimNode
{
    Memory &memory;
    U letter;
    U sr1;
    U sr2;
    Do *E;
    Do *F;
    Do *G;
    int Efn;
    int Ffn;
    int Gfn;
    char *now;
    size_t remain;
    size_t chunk_id;

    SimNode(Memory &memory_, U letter_, U sr1p, U sr2p, int Efn_, int Ffn_, int Gfn_)
        : memory(memory_)
    {
        letter = letter_;
        edge.Bro.PreRange(sr1p, sr2p, letter);
        sr1 = sr1p;
        sr2 = sr2p;

        now = memory.now;
        remain = memory.remain;
        chunk_id = memory.chunk_id;
        Efn = Efn_;
        Ffn = Ffn_;
        Gfn = Gfn_;
        E = memory.heap_alloc<Do>(Efn);
        F = memory.heap_alloc<Do>(Ffn);
        G = memory.heap_alloc<Do>(Gfn);
    }

    ~SimNode()
    {
        memory.now = now;
        memory.remain = remain;
        memory.chunk_id = chunk_id;
    }
};

struct Dot // cannot use constructor because on Memory
{
    int n; // n determines the type of Dot
    size_t s;
    int w;
    double val;
    Dot **sources; // apply on memory
    int s_sz;      // s_sz<0 means visited
};

struct Ptr
{
    Dot Abar;
    Dot *B;
    Dot *C;
    Dot *D;
    Dot **A;
    Dot **E;
    Dot **F;
    Dot **G;
    Dot **F0;
    Dot **G0;
    Dot **D0;
    Dot **DX;
    std::set<std::pair<size_t,int>> Cdelta;
    std::map<int,std::set<std::pair<size_t,int>>> Adelta;
};

template <typename T>
struct SNC
{
    double E;
    double F0;
    double G0;
    double hatG;
    T sr1;
    T sr2;
    SNC *itp;
    SNC *itcs[5]; // NULL means children are not added to the tree yet
};

template <typename U>
struct Align : Memory
{
    const static int ww = 3;

    Graph<U> &graph;
    std::stack<SimNode<U>> simnodes;
    std::map<std::string, Ptr> name_ptr;
    std::string O;

    Align(Graph<U> &graph_, size_t chunk_sz_) : graph(graph_)
    {
        Initial(chunk_sz_);
        for (auto &node : graph.nodes)
        {
            Ptr &ptr = name_ptr[node.first];
            if (graph.roots.find(node.first) == graph.roots.end())
                ptr.Abar.val = -inf;
            else
                ptr.Abar.val = 0;
            ptr.Abar.n = -3;
            ptr.Abar.s = 0;
            ptr.Abar.w = 0;
            ptr.Abar.s_sz = 0;
            ptr.Abar.sources = NULL;
        }
    }

    int base2int(char c)
    {
        return (c + 6) >> 1 & 7;
    }

    void Mix()
    {
        clear();
        for (auto &node : graph.nodes)
        {
            auto &ptr = name_ptr[node.first];
            alloc_initial(ptr.A, -1, node.second.scc_ptr->node_names.size() - 1);
            alloc_initial(ptr.B, -2, 0);
        }
        for (auto &scc : graph.sccs)
        {
            if (scc.node_names.size() != 0)
            {
                for (auto &edge_name : scc.local_names)
                {
                    auto &ptr = name_ptr[edge_name];
                    auto &edge = graph.locals[edge_name];
                    alloc_initial(ptr.D0, edge.n, scc.node_names.size() - 1);
                    alloc_initial(ptr.DX, edge.n, scc.node_names.size() - 1);
                    alloc_initial(ptr.D, edge.n, 0);
                    alloc_initial(ptr.E, edge.n, edge.seq.size());
                    alloc_initial(ptr.F0, edge.n, edge.seq.size());
                    alloc_initial(ptr.G0, edge.n, edge.seq.size());
                    alloc_initial(ptr.G, edge.n, edge.seq.size());
                }
                for (auto &edge_name : scc.global_names)
                {
                    auto &ptr=name_ptr[edge_name];
                    auto &edge = graph.globals[edge_name];
                    alloc_initial(ptr.C, edge.n, 0);
                    alloc_initial(ptr.D0, edge.n, scc.node_names.size() - 1);
                }
                for (int w=0; w<=O.size(); ++w)
                {
                    for (auto &edge_name : scc.local_names)
                        CircuitIteration(w, edge_name);
                    for (auto &edge_name : scc.global_names)
                        CircuitIterationGlobal(w, edge_name);
                    for (auto &node_name : scc.node_names)
                    {
                        auto &ptr = name_ptr[node_name];
                        auto &node = nodes[node_name];
                        if (w == 0)
                            source_max(ptr.B[w],{},{});
                        else
                            source_max(ptr.B[w],{&ptr.B[w-1],&ptr.A[scc.node_names.size()-1][w-1]},{ptr.B[w-1].val+node.ue,ptr.A[scc.node_names.size()-1][w-1]+node.ve});
                        std::vector<Dot *> adrs{&ptr.B[w]};
                        std::vector<double> vals{ptr.B[w].val};
                        for (auto &edge_name : node.in_edges)
                            
                        source_max(ptr.A[0][w],{});
                    }
                }
            }
            else
            {

            }
        }
        for (auto &edge_name : graph.edge_names)
        {
            auto &ptr = name_ptr[edge_name];
            auto &edge = graph.globals.find(edge_name) == graph.globals.end() ? graph.locals[edge_name] : graph.globals[edge_name];
            int scc_sz = edge.scc_ptr->node_names.size();
            alloc_initial(ptr.D0, edge.n, scc_sz - 1);
            if (graph.globals.find(edge_name) != graph.globals.end())
                continue;
            alloc_initial(ptr.DX, edge.n, scc_sz - 1);
            if (scc_sz == 0)
            {
                alloc_initial(ptr.E, edge.n, edge.seq.size());
                alloc_initial(ptr.F, edge.n, edge.seq.size());
                alloc_initial(ptr.G, edge.n, edge.seq.size());
            }
            else
            {
                alloc_initial(ptr.D, edge.n, 0);
                alloc_initial(ptr.E, edge.n, edge.seq.size());
                alloc_initial(ptr.F0, edge.n, edge.seq.size());
                alloc_initial(ptr.G0, edge.n, edge.seq.size());
                alloc_initial(ptr.G, edge.n, edge.seq.size());
            }
        }
        
    }

    void CrossIteration(std::string &edge_name)
    {
        EdgeLocal &edge = graph.locals[edge_name];
        Dot **A = name_prt[edge.tail].A;
        Dot **&E = name_prt[edge_name].E;
        Dot **&F = name_prt[edge_name].F;
        Dot **&G = name_prt[edge_name].G;
        int tail_scc_sz = graph.nodes[edge.tail].scc_ptr->node_names.size();
        for (size_t s = 0; s <= edge.seq.size(); ++s)
            for (size_t w = 0; w <= O.size(); ++w)
            {
                if (w == 0)
                    source_max(E[s][w], {}, {});
                else
                    source_max(E[s][w], {&E[s][w - 1], &G[s][w - 1]}, {E[s][w - 1].val + edge.ue, G[s][w - 1].val + edge.ve});
                if (s == 0)
                    source_max(F[s][w], {}, {});
                else
                    source_max(F[s][w], {&F[s - 1][w], &G[s - 1][w]}, {F[s - 1][w].val + edge.uf, G[s - 1][w].val + edge.vf});
                if (s > 0 && w > 0)
                    source_max(G[s][w], {&F[s][w], &E[s][w], &A[tail_scc_sz - 1][w]}, {F[s][w].val, E[s][w].val, A[tail_scc_sz - 1][w].val + (s == 0 ? 0 : edge.vfp + (s - 1) * edge.ufp) + edge.T});
                else
                    source_max(G[s][w], {&F[s][w], &E[s][w], &G[s - 1][w - 1], &A[tail_scc_sz - 1][w]}, {F[s][w].val, E[s][w].val, G[s - 1][w - 1].val + edge.gamma[base2int(edge.seq[s])][base2int(O[w])], A[tail_scc_sz - 1][w].val + (s == 0 ? 0 : edge.vfp + (s - 1) * edge.ufp) + edge.T});
            }
    }

    void CrossIterationGlobal(int g, std::string &O)
    {
        EdgeGlobal<U> &edge = graph.globals[g];
        EdgeGlobalCopy<U> &edgecopy = this->globalcopys[g];
        Dot<U> *tildeANT = this->nodecopys[edge.tail].tildeA[graph.nodes[edge.tail].tAsz - 1];
        if (tildeANT[0].val == -inf)
            return;
        Node<U> &NH = graph.nodes[edge.head];
        double tve = NH.tve, tue = NH.tue;
        BroWheel<U> &Bro = edge.Bro;
        double ve = edge.ve, ue = edge.ue, T = edge.T;
        std::vector<std::vector<U>> &boxes = edgecopy.boxes;

        Dot<U> *A = edgecopy.A[0], *B = edgecopy.B;
        for (size_t w = 0; w <= O.size(); ++w)
        {
            A[w].val = -inf;
            boxes[w].resize(3);
            boxes[w][0] = boxes[w][2] = 0;
            boxes[w][1] = Bro.size() - 1;
        }
        simnodes.emplace();
        simnodes.top().sr1 = 0;
        simnodes.top().sr2 = Bro.size() - 1;

        triestates.emplace();
        TrieState<U> &ts = triestates.top();
        PushTrieState(ts, O.size() + 1, 0, O.size() + 1);
        ts.E[0].w = ts.G[0].w = 0;
        ts.E[0].val = -inf;
        ts.G[0].val = tildeANT[0].val + T;
        A[0].val = ts.G[0].val;
        for (size_t g = 1; g <= O.size(); ++g)
        {
            ts.G[g].w = ts.E[g].w = g;
            ts.E[g].val = std::max(ts.E[g - 1].val + ue, ts.G[g - 1].val + ve);
            ts.G[g].val = std::max(ts.E[g].val, tildeANT[g].val + T);
            A[g].val = ts.G[g].val;
        }
        simnodes.top().visit = true;
        PushPre(edge);

        while (!simnodes.empty())
        {
            if (simnodes.top().visit)
            {
                TrieState<U> &ts = triestates.top();
                this->now = ts.now;
                this->remain = ts.remain;
                this->chunk_id = ts.chunk_id;
                triestates.pop();
                simnodes.pop();
            }
            else
            {
                TrieState<U> &tsp = triestates.top();
                triestates.emplace();
                TrieState<U> &ts = triestates.top();
                PushTrieState(ts, O.size() - tsp.G[0].w + 1, tsp.Gfn, O.size() - tsp.G[0].w + 1);
                CrossBodyGlobal(tsp, ts, edge, edgecopy, O, simnodes.top().letter, simnodes.top().sr1, simnodes.top().sr2, triestates.size() - 1);
                simnodes.top().visit = true;
                if (ts.Gfn > 0)
                {
                    if (simnodes.top().sr1 < simnodes.top().sr2)
                        PushPre(edge);
                    else
                    {
                        U os = edge.Bro.SimSuffix(simnodes.top().sr1);
                        if (--os >= 0)
                        {
                            TrieState<U> tssp, tss;
                            PushTrieState(tssp, O.size() + 1, O.size() + 1, O.size() + 1);
                            PushTrieState(tss, O.size() - ts.G[0].w + 1, ts.Gfn, O.size() - ts.G[0].w + 1);
                            U letter = edge.G(os);
                            U sr = simnodes.top().sr1;
                            sr = Bro.C[letter - 1] + Bro.GetOcc(sr, letter);
                            U seq_sz = triestates.size();
                            CrossBodyGlobal(ts, tss, edge, edgecopy, O, letter, sr, sr, seq_sz);
                            while (tss.Gfn > 0 && --os >= 0)
                            {
                                std::swap(tssp, tss);
                                this->now = tss.now, this->remain = tss.remain, this->chunk_id = tss.chunk_id;
                                tss.Efn = tss.Gfn = O.size() - tssp.G[0].w + 1;
                                tss.Ffn = tssp.Gfn;
                                tss.E = heap_alloc<Do>(tss.Efn);
                                tss.F = heap_alloc<Do>(tss.Ffn);
                                tss.G = heap_alloc<Do>(tss.Gfn);
                                letter = edge.G(os);
                                sr = Bro.C[letter - 1] + Bro.GetOcc(sr, letter);
                                CrossBodyGlobal(tssp, tss, edge, edgecopy, O, letter, sr, sr, ++seq_sz);
                            }
                        }
                    }
                }
            }
        }
        B[0].val = -inf;
        A[0].s_sz = 0;
        for (size_t w = 1; w <= O.size(); ++w)
        {
            source_max(B[w], {&B[w - 1], &A[w - 1]}, {B[w - 1].val + tue, A[w - 1].val + tve});
            if (B[w].val >= A[w].val)
            {
                if (B[w].val > A[w].val)
                    boxes[w].clear();
                source_max(A[w], {&B[w]}, {B[w].val});
            }
            else
                A[w].s_sz = 0;
        }
        ++this->nodecopys[edge.head].cid_now;
        update_cross(this->nodecopys[edge.head], O.size());
    }

    void CrossBodyGlobal(TrieState<U> &tsp, TrieState<U> &ts, EdgeGlobal<U> &edge, EdgeGlobalCopy<U> &edgecopy, std::string &O, U letter, U sr1, U sr2, U seq_sz)
    {
        Dot<U> *tildeANT = this->nodecopys[edge.tail].tildeA[graph.nodes[edge.tail].tAsz - 1];
        for (int g = 0, f = 0; g < ts.Ffn; ++g)
        {
            int w = ts.F[g].w = tsp.G[g].w;
            ts.F[g].val = tsp.G[g].val + edgecopy.vf[w];
            if (f < tsp.Ffn && w == tsp.F[f].w)
            {
                ts.F[g].val = std::max(ts.F[g].val, tsp.F[f].val + edgecopy.uf[w]);
                ++f;
            }
        }
        ts.E[0].w = ts.G[0].w = O.size() - ts.Gfn + 1;
        ts.E[0].val = -inf;
        if (ts.F[0].val < tildeANT[ts.G[0].w].val + edge.T)
            ts.F[0].val = -inf;
        ts.G[0].val = ts.F[0].val;
        for (int g = 1, f = 1; g < ts.Gfn; ++g)
        {
            int w = ts.E[g].w = ts.G[g].w = ts.G[g - 1].w + 1;
            ts.E[g].val = std::max(ts.E[g - 1].val + edge.ue, ts.G[g - 1].val + edge.ve);
            if (ts.E[g].val + edge.ue < tildeANT[w].val + edge.T + edge.ve)
                ts.E[g].val = -inf;
            ts.G[g].val = ts.E[g].val;
            if (tsp.G[f - 1].w == w - 1)
                ts.G[g].val = std::max(ts.G[g].val, tsp.G[f - 1].val + edge.gamma[letter][int(O[w - 1])]);
            if (f < ts.Ffn && ts.F[f].w == w)
            {
                ts.G[g].val = std::max(ts.G[g].val, ts.F[f].val);
                ++f;
            }
            if (ts.G[g].val < tildeANT[w].val + edge.T)
            {
                ts.G[g].val = -inf;
                if (ts.F[f - 1].w == w)
                    ts.F[f - 1].val = -inf;
            }
        }
        {
            int *Xfns[2] = {&ts.Ffn, &ts.Gfn};
            int i = 0;
            for (Do *X : {ts.F, ts.G})
            {
                int XfnNew = 0;
                for (int x = 0; x < *Xfns[i]; ++x)
                    if (X[x].val != -inf)
                    {
                        if (x != XfnNew)
                            X[XfnNew] = X[x];
                        ++XfnNew;
                    }
                *Xfns[i] = XfnNew;
                ++i;
            }
        }
        for (int g = 0; g < ts.Gfn; ++g)
        {
            int w = ts.G[g].w;
            if (ts.G[g].val >= edgecopy.A[0][w].val)
            {
                if (ts.G[g].val > edgecopy.A[0][w].val)
                {
                    edgecopy.A[0][w].val = ts.G[g].val;
                    edgecopy.boxes[w].clear();
                }
                edgecopy.boxes[w].push_back(sr1);
                edgecopy.boxes[w].push_back(sr2);
                edgecopy.boxes[w].push_back(seq_sz);
            }
        }
    }

    void CircuitIteration(int w, std::string &edge_name)
    {
        EdgeLocal &edge = graph.locals[edge_name];
        Dot **A = name_ptr[edge.tail].A;
        Ptr &ptr = name_ptr[edge_name];
        int tail_scc_sz = edge.scc_ptr->node_names.size();
        if (w > 0)
    }

    template <typename Y>
    void CircuitInitial(int w, Dot<U> **C, Dot<U> **E, Dot<U> **F, Dot<U> **G, Dot<U> **H, double ve, double ue, Y &edge)
    {
        Dot<U> *tildeANT = this->nodecopys[edge.tail].tildeA[graph.nodes[edge.tail].tAsz - 1];
        C[0][w - 1].val = -inf;
        source_max(H[0][w - 1], {&tildeANT[w - 1]}, {tildeANT[w - 1].val + edge.T});
        source_max(E[0][w], {&E[0][w - 1], &H[0][w - 1]}, {E[0][w - 1].val + ue, H[0][w - 1].val + ve});
        F[0][w].val = -inf;
        source_max(G[0][w], {&E[0][w]}, {E[0][w].val});
    }

    template <typename Y, typename Z>
    void CircuitBody(std::string &seq, std::string &O, int w, int s, Dot<U> **C, Dot<U> **E, Dot<U> **F, Dot<U> **G, Dot<U> **H, double ve, double ue, Y &edge, Z &edgecopy, double gfps)
    {
        Dot<U> *tildeANT = this->nodecopys[edge.tail].tildeA[graph.nodes[edge.tail].tAsz - 1];
        source_max(C[s][w - 1], {&C[s - 1][w - 1], &H[s - 1][w - 1]}, {C[s - 1][w - 1].val + edgecopy.uf[w - 1], H[s - 1][w - 1].val + edgecopy.vf[w - 1]});
        source_max(H[s][w - 1], {&C[s][w - 1], &tildeANT[w - 1]}, {C[s][w - 1].val, tildeANT[w - 1].val + gfps + edge.T});
        source_max(E[s][w], {&G[s][w - 1], &E[s][w - 1], &H[s][w - 1]}, {G[s][w - 1].val + ve, E[s][w - 1].val + ue, H[s][w - 1].val + ve});
        source_max(F[s][w], {&G[s - 1][w], &F[s - 1][w]}, {G[s - 1][w].val + edgecopy.vf[w], F[s - 1][w].val + edgecopy.uf[w]});
        double gm = edge.gamma[int(seq[s - 1])][int(O[w - 1])];
        source_max(G[s][w], {&E[s][w], &F[s][w], &G[s - 1][w - 1], &H[s - 1][w - 1]}, {E[s][w].val, F[s][w].val, G[s - 1][w - 1].val + gm, H[s - 1][w - 1].val + gm});
    }

    void CircuitIterationEdgeGlobal(int w, int g, std::string &O)
    {
        EdgeGlobal<U> &edge = graph.globals[g];
        EdgeGlobalCopy<U> &edgecopy = this->globalcopys[g];
        NodeCopy<U> &NTCP = this->nodecopys[edge.tail];
        NodeCopy<U> &NHCP = this->nodecopys[edge.head];
        std::forward_list<SNC<U> *> &sncs = edgecopy.sncs;
        Dot<U> *B = edgecopy.B, **A = edgecopy.A;
        BroWheel<U> &Bro = edge.Bro;
        double T = edge.T, ue = edge.ue, ve = edge.ve;
        std::vector<double> &uf = edgecopy.uf, &vf = edgecopy.vf;
        std::vector<std::vector<U>> &boxes = edgecopy.boxes;
        if (w == 0)
        {
            sncs.emplace_front();
            sncs.front() = memtree.heap_alloc<SNC<U>>(1);
            sncs.front()->E = sncs.front()->F = sncs.front()->G = -inf;
            sncs.front()->s = 0;
            sncs.front()->sr1 = 0;
            sncs.front()->sr2 = Bro.size() - 1;
            sncs.front()->itcs = NULL;
            B[0].val = -inf;
            A[0][w].val = -inf;
        }
        else
        {
            double H = NTCP.tildeA[graph.nodes[edge.tail].tAsz - 1][w - 1].val + T;
            auto lit = sncs.begin();
            SNC<U> *it = *lit;
            it->Gp = it->G;
            it->F = -inf;
            A[0][w].val = it->G = it->E = std::max(it->E + ue, H + ve);
            boxes[w].clear();
            boxes[w].push_back(it->sr1);
            boxes[w].push_back(it->sr2);
            boxes[w].push_back(it->s);
            InsertSNC(lit, sncs, edge);
            while (lit != sncs.end())
            {
                it = *lit;
                it->Gp = it->G;
                it->E = std::max(it->G + ve, it->E + ue);
                it->F = std::max(it->itp->G + vf[w], it->itp->F + uf[w]);
                double gm = edge.gamma[it->letter][int(O[w - 1])];
                it->G = std::max(std::max(it->E, it->F), it->itp->Gp + gm);
                if (it->s == 1)
                    it->G = std::max(it->G, H + gm);
                if (it->G < NTCP.tildeA[0][w].val + T)
                    it->G = it->F = -inf;
                if (it->E + ue < NTCP.tildeA[0][w].val + T + ve)
                    it->E = -inf;
                if (it->G >= A[0][w].val)
                {
                    if (it->G > A[0][w].val)
                    {
                        A[0][w].val = it->G;
                        boxes[w].clear();
                    }
                    boxes[w].push_back(it->sr1);
                    boxes[w].push_back(it->sr2);
                    boxes[w].push_back(it->s);
                }
                InsertSNC(lit, sncs, edge);
            }
            if (size_t(w) == O.size())
                sncs.clear();
            else
            {
                auto lit = sncs.begin(), nlit = std::next(lit);
                while (nlit != sncs.end())
                {
                    SNC<U> *nit = *nlit;
                    if (nit->G == -inf && nit->E == -inf && nit->itp->G == -inf)
                    {
                        nit->Gp = -inf;
                        sncs.erase_after(lit);
                    }
                    else
                        ++lit;
                    nlit = std::next(lit);
                }
            }
            if (B[w].val >= A[0][w].val)
            {
                if (B[w].val > A[0][w].val)
                    boxes[w].clear();
                source_max(A[0][w], {&B[w]}, {B[w].val});
            }
            else
                A[0][w].s_sz = 0;
            if (A[0][w].val >= NHCP.tildeA[0][w].val)
            {
                adrs.clear();
                vals.clear();
                for (int j = 0; j < NHCP.tildeA[0][w].s_sz; ++j)
                {
                    adrs.push_back(NHCP.tildeA[0][w].sources[j]);
                    vals.push_back(NHCP.tildeA[0][w].val);
                }
                adrs.push_back(&A[0][w]);
                vals.push_back(A[0][w].val);
                source_max(NHCP.tildeA[0][w], adrs.begin(), vals.begin(), vals.end());
            }
        }
    }

    void InsertSNC(typename std::forward_list<SNC<U> *>::iterator &lit, std::forward_list<SNC<U> *> &sncs, EdgeGlobal<U> &edge)
    {
        SNC<U> *it = *lit;
        if (it->G != -inf && it->Gp == -inf)
        {
            if (it->itcs == NULL)
            {
                if (it->sr1 == it->sr2)
                {
                    if (it->itp->sr1 < it->itp->sr2)
                        it->itcs_sz = edge.Bro.SimSuffix(it->sr1);
                    if (it->itcs_sz > 0)
                    {
                        SNC<U> *nit = memtree.heap_alloc<SNC<U>>(1);
                        nit->itcs_sz = it->itcs_sz - 1;
                        nit->s = it->s + 1;
                        nit->letter = edge.G(nit->itcs_sz);
                        nit->sr1 = nit->sr2 = edge.Bro.C[nit->letter - 1] + edge.Bro.GetOcc(it->sr1, nit->letter);
                        nit->E = nit->F = nit->G = -inf;
                        nit->itp = it;
                        nit->itcs = NULL;
                        it->itcs = memtree.heap_alloc<SNC<U> *>(1);
                        it->itcs[0] = nit;
                    }
                    else
                        it->itcs = (SNC<U> **)(memtree.now);
                }
                else
                {
                    U left = it->sr2 - it->sr1 + 1;
                    std::queue<SNC<U> *> itcs;
                    for (U letter = 2; letter < 7; ++letter)
                    {
                        U sr1 = it->sr1, sr2 = it->sr2;
                        edge.Bro.PreRange(sr1, sr2, letter);
                        if (sr1 <= sr2)
                        {
                            SNC<U> *nit = memtree.heap_alloc<SNC<U>>(1);
                            itcs.push(nit);
                            nit->sr1 = sr1;
                            nit->sr2 = sr2;
                            nit->s = it->s + 1;
                            nit->letter = letter;
                            nit->E = nit->F = nit->G = -inf;
                            nit->itp = it;
                            nit->itcs = NULL;
                            left -= nit->sr2 - nit->sr1 + 1;
                            if (left == 0)
                                break;
                        }
                    }
                    it->itcs_sz = itcs.size();
                    it->itcs = memtree.heap_alloc<SNC<U> *>(it->itcs_sz);
                    for (int i = 0; i < it->itcs_sz; ++i)
                    {
                        it->itcs[i] = itcs.front();
                        itcs.pop();
                    }
                }
            }
            if (it->itcs_sz != 0)
            {
                size_t sz;
                if (it->sr1 == it->sr2)
                    sz = 1;
                else
                    sz = it->itcs_sz;
                for (size_t i = 0; i < sz; ++i)
                {
                    SNC<U> *itc = it->itcs[i];
                    if (itc->E == -inf && itc->G == -inf)
                        sncs.emplace_after(lit, itc);
                }
            }
        }
        ++lit;
    }

    void update_cross(NodeCopy<U> &nodecopy, int W)
    {
        if (size_t(nodecopy.cid_now) == nodecopy.ciptr.size())
            for (int w = 0; w <= W; ++w)
            {
                adrs.clear();
                vals.clear();
                if (nodecopy.barA)
                {
                    adrs.push_back(&nodecopy.barA[w]);
                    vals.push_back(nodecopy.barA[w].val);
                }
                for (Dot<U> **ptr : nodecopy.ciptr)
                {
                    adrs.push_back(&ptr[0][w]);
                    vals.push_back(ptr[0][w].val);
                }
                source_max(nodecopy.tildeA[0][w], adrs.begin(), vals.begin(), vals.end());
            }
    }

    void alloc_initial(Dot **ptr, int n, int S)
    {
        ptr = heap_alloc<Dot *>(S + 1);
        for (int s = 0; s <= S; ++s)
            alloc_initial(ptr[s], int n, int s);
    }

    void alloc_initial(Dot *ptr, int n, int s)
    {
        ptr = heap_alloc<Dot>(O.size() + 1);
        for (int w = 0; w <= O.size(); ++w)
        {
            ptr[w].n = n;
            ptr[w].s = s;
            ptr[w].w = w;
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
        for (Itv itv = itv_b, int j = 0; itv != itv_e; ++itv, ++ita_b, ++j)
            if (*itv == dot.val)
                dot.sources[j] = *ita_b;
    }
};

#endif
