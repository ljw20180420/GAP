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

struct Align : Graph
{
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

        

        struct Do
        {
            QUERYSIZE w;
            SCORETYPE E;
        };

        std::unique_ptr<Do []> E;
    };

    std::mutex &mtx;
    std::ifstream &fin;
    std::vector<NUCTYPE> O;
    std::string Oname;
    std::ofstream fout;
    IDTYPE max_id = 0;
    QUERYSIZE Omax;
    std::unique_ptr<IDTYPE []> ids;
    std::unique_ptr<BITSTYPE []> bits;

    CrossGlobalData crossglobaldata;

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

    std::unique_ptr<SCORETYPE []> croS;
    SCORETYPE *Ethres, *croE;
    // std::unique_ptr<SCORETYPE *[]> croF, croG;
    // std::unique_ptr<QUERYSIZE []> croQ;
    // std::unique_ptr<QUERYSIZE *[]> crop;

    SIZETYPE trn, tnn;

    Align(std::mutex &mtx_, std::ifstream &fin_, boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short, std::map<std::string, RankVec> &file2rankvec, std::string mgfile, QUERYSIZE Omax_)
        : mtx(mtx_), fin(fin_), Graph(vm, file2short, file2rankvec), fout(mgfile, std::ifstream::binary), Omax(Omax_)
    {
        trn = 0;
        tnn = 1 + nodes.size(); // Q and Abars
        for (Node &node : nodes)
        {
            trn += node.scc_sz + 1;
            tnn += (Omax + 1) * (node.scc_sz + 1);
        }
        for (EdgeLocalCross &edge : local_crosses)
        {
            SHORTSIZE ref_sz = file2short[edge.name].second;
            trn += 3 * (ref_sz + 1);
            tnn += (Omax + 1) * 3 * (ref_sz + 1);
        }
        for (EdgeLocalCircuit &edge : local_circuits)
        {
            SHORTSIZE ref_sz = file2short[edge.name].second;
            trn += 1 + 4 * (ref_sz + 1) + 2 * (edge.tail->scc_sz - 1);
            tnn += (Omax + 1) * (1 + 4 * (ref_sz + 1) + 2 * (edge.tail->scc_sz - 1));
        }
        for (EdgeGlobalCircuit &edge : global_circuits)
        {
            trn += edge.tail->scc_sz - 1;
            tnn += (Omax + 1) * (edge.tail->scc_sz - 1);
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
        for (SIZETYPE i=0; i<nodes.size(); ++i)
        {
            nodes[i].pAbarval = --rpval;
            *(fps++) = 0;
            *(fpn++) = Dot::NODEIDX2DOTTYPEB(nodes[i].n);
            for (SIZETYPE j = 0; j < nodes[i].scc_sz; ++j, ++fps, ++fpn)
            {
                *fps = j;
                *fpn = Dot::NODEIDX2DOTTYPEA(nodes[i].n);
            }
            nodes[i].Bvals = fpval;
            fpval += Omax + 1;
            nodes[i].Avals.reset(new SCORETYPE *[nodes[i].scc_sz]);
            for (SIZETYPE s = 0; s < nodes[i].scc_sz; ++s, fpval += Omax + 1)
                nodes[i].Avals[s] = fpval;
            nodes[i].AdeltaDot.reset(new std::deque<SCORETYPE *>[Omax + 1]);
            nodes[i].AdeltaGlobal.reset(new std::deque<Node::GlobalSuffix>[Omax + 1]);
        }

        for (EdgeLocalCross &edge : local_crosses)
        {
            SHORTSIZE ref_sz = file2short[edge.name].second;
            for (SIZETYPE i = 0; i < 3; ++i)
                for (SIZETYPE j = 0; j <= ref_sz; ++j, ++fps, ++fpn)
                {
                    *fps = j;
                    *fpn = edge.n;
                }
            SIZETYPE si = 0;
            for (std::unique_ptr<SCORETYPE *[]> *pYvalss : {&edge.Evals, &edge.Fvals, &edge.Gvals})
            {
                pYvalss->reset(new SCORETYPE *[ref_sz + 1]);
                for (SIZETYPE s = 0; s <= ref_sz; ++s, fpval += Omax + 1)
                    (*pYvalss)[s] = fpval;
            }
        }

        for (EdgeLocalCircuit &edge : local_circuits)
        {
            SHORTSIZE ref_sz = file2short[edge.name].second;
            *(fps++) = 0;
            *(fpn++) = edge.n;
            SIZETYPE Ss[6] = {ref_sz+1, ref_sz+1, ref_sz+1, ref_sz+1, edge.tail->scc_sz-1, edge.tail->scc_sz-1};
            for (SIZETYPE i = 0; i < 6; ++i)
                for (SIZETYPE j = 0; j < Ss[i]; ++j, ++fps, ++fpn)
                {
                    *fps = j;
                    *fpn = edge.n;
                }
            edge.Dvals = fpval;
            fpval += Omax + 1;
            SIZETYPE si = 0;
            for (std::unique_ptr<SCORETYPE *[]> *pYvalss : {&edge.Evals, &edge.F0vals, &edge.G0vals, &edge.Gvals, &edge.D0vals, &edge.DXvals})
            {
                pYvalss->reset(new SCORETYPE *[Ss[si]]);
                for (SIZETYPE s = 0; s < Ss[si]; ++s, fpval += Omax + 1)
                    (*pYvalss)[s] = fpval;
                ++si;
            }
        }

        for (EdgeGlobalCircuit &edge : global_circuits)
        {
            for (SIZETYPE j = 0; j < edge.tail->scc_sz-1; ++j, ++fps, ++fpn)
            {
                *fps = j;
                *fpn = edge.n;
            }
            edge.D0vals.reset(new SCORETYPE *[edge.tail->scc_sz-1]);
            for (SIZETYPE s = 0; s < edge.tail->scc_sz-1; ++s, fpval += Omax + 1)
                edge.D0vals[s] = fpval;
        }

        crossglobaldata.E.reset(new CrossGlobalData::Do[Omax + 1]);

        // croS.reset(new SCORETYPE[(2 + 4 * Omax) * (Omax + 1)]);
        croS.reset(new SCORETYPE[2 * (Omax + 1)]);
        Ethres = croS.get();
        croE = Ethres + Omax + 1;
        // SCORETYPE *ptr = croE + Omax + 1;
        // for (std::unique_ptr<SCORETYPE *[]> *pcroX : {&croF, &croG})
        // {
        //     pcroX->reset(new SCORETYPE *[2 * Omax]);
        //     for(SIZETYPE i = 0; i < 2 * Omax; ++i, ptr += Omax + 1)
        //         (*pcroX)[i] = ptr;
        // }
        // croQ.reset(new QUERYSIZE[2 * Omax * (Omax + 2)]);
        // ptr = croQ.get();
        // crop.reset(new QUERYSIZE *[2 * Omax]);
        // for(SIZETYPE i = 0; i < 2 * Omax; ++i, ptr += Omax + 2)
        //     crop[i] = ptr;
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
        for (Node &node : nodes)
        {
            for (SIZETYPE w = 0; w <= O.size(); ++w)
                node.Avals[0][w] = -inf;
            if (node.is_root)
            {
                *node.pAbarval = 0;
                node.Avals[0][0] = 0;
                node.AdeltaDot[0].emplace_back(node.pAbarval);
            }
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

                    SIZETYPE M = &node.Bvals[w] - vals.get();
                    bits[M] = 0;
                    if (w == 0)
                        node.Bvals[w] = -inf;
                    else
                    {
                        node.Bvals[w] = std::max(node.Bvals[w - 1] + node.ue, node.Avals[node.scc_sz - 1][w - 1] + node.ve);
                        if (node.Bvals[w] == node.Bvals[w - 1] + node.ue)
                            bits[M] += 1;
                        if (node.Bvals[w] == node.Avals[node.scc_sz - 1][w - 1] + node.ve)
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
                        SIZETYPE M = &(circuit->DXvals[l - 1][w]) - vals.get();
                        bits[M] = 0;
                        circuit->DXvals[l - 1][w] = std::max(circuit->tail->Avals[l - 1][w] + circuit->gfpT[ref_sz], circuit->D0vals[l - 1][w] + circuit->gf[ref_sz]);
                        if (circuit->DXvals[l - 1][w] == circuit->tail->Avals[l - 1][w] + circuit->gfpT[ref_sz])
                            bits[M] += 1;
                        if (circuit->DXvals[l - 1][w] == circuit->D0vals[l - 1][w] + circuit->gf[ref_sz])
                            bits[M] += 2;
                    }
                    for (Node &node : nodes)
                    {
                        if (node.scc_id != scc_id)
                            continue;
                        
                        SIZETYPE M = &node.Avals[l][w] - vals.get();
                        node.Avals[l][w] = node.Avals[0][w];
                        bits[M] = 1;
                        SIZETYPE ei = 1;
                        for (Edge *edge : edges)
                        {
                            if (edge->head != &node || edge->tail->scc_id != scc_id)
                                continue;
                            ei *= 2;
                            if (ei > std::numeric_limits<BITSTYPE>::max())
                            {
                                std::cerr << "the upper limits of BITSTYPE reached by node" << node.n << ", A[" << l << "][" << w << "]\n";
                                break;
                            }
                            SCORETYPE Dscore;
                            if (file2short.count(edge->name))
                            {
                                EdgeLocalCircuit *circuit = (EdgeLocalCircuit *)edge;
                                Dscore = std::max(circuit->D0vals[l - 1][w] + circuit->gfm, circuit->DXvals[l - 1][w]);
                            }
                            else
                            {
                                EdgeGlobalCircuit *circuit = (EdgeGlobalCircuit *)edge;
                                Dscore = circuit->D0vals[l - 1][w];
                            }
                            if (Dscore >= node.Avals[l][w])
                            {
                                if (Dscore > node.Avals[l][w])
                                {
                                    node.Avals[l][w] = Dscore;
                                    bits[M] = 0;
                                }
                                bits[M] += ei;
                            }
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

        SIZETYPE M = pQval - vals.get();
        bits[M] = 0;
        *pQval = -inf;
        SIZETYPE ei = 1;
        for (Node &node : nodes)
        {
            if (!node.is_target)
                continue;
            if (ei > std::numeric_limits<BITSTYPE>::max())
            {
                std::cerr << "the upper limits of BITSTYPE reached by Q (the graph has too many targets)\n";
                break;
            }
            SCORETYPE score = node.Avals[node.scc_sz - 1][O.size()];
            if (score >= *pQval)
            {
                if (score > *pQval)
                {
                    *pQval = score;
                    bits[M] = 0;
                }
                bits[M] += ei;
            }
            ei *= 2;
        }
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
                SIZETYPE M = &Evals[s][w] - vals.get();
                bits[M] = 0;
                if (w == 0)
                    vals[M] = -inf;
                else
                {
                    vals[M] = std::max(Evals[s][w - 1] + edge->ue, Gvals[s][w - 1] + edge->ve);
                    if (vals[M] == Evals[s][w - 1] + edge->ue)
                        bits[M] += 1;
                    if (vals[M] == Gvals[s][w - 1] + edge->ve)
                        bits[M] += 2;
                }

                M = &Fvals[s][w] - vals.get();
                bits[M] = 0;
                if (s == 0)
                    vals[M] = -inf;
                else
                {
                    vals[M] = std::max(Fvals[s - 1][w] + edge->uf, Gvals[s - 1][w] + edge->vf);
                    if (vals[M] == Fvals[s - 1][w] + edge->uf)
                        bits[M] += 1;
                    if (vals[M] == Gvals[s - 1][w] + edge->vf)
                        bits[M] += 2;
                }
                
                SCORETYPE scores[4] = {Evals[s][w], Fvals[s][w], Avals[w] + edge->gfpT[s], -inf};
                if (s > 0 && w > 0)
                    scores[3] = Gvals[s - 1][w - 1] + edge->gamma[ref[s - 1]][O[w - 1]];
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


                if (Gvals[s][w] >= edge->head->Avals[0][w])
                {
                    if (Gvals[s][w] > edge->head->Avals[0][w])
                    {
                        edge->head->Avals[0][w] = Gvals[s][w];
                        edge->head->AdeltaDot[w].clear();
                        edge->head->AdeltaGlobal[w].clear();
                    }
                    edge->head->AdeltaDot[w].emplace_back(&Gvals[s][w]);
                }
            }
    }

    void CrossIterationGlobal(Edge *edge)
    {
        SCORETYPE *Avals = edge->tail->Avals[edge->tail->scc_sz - 1];
        std::unique_ptr<CrossGlobalData::Do[]> &E = crossglobaldata.E;
        std::deque<CrossGlobalData::Doo> &FG = crossglobaldata.FG;
        RankVec &rankvec = file2rankvec[edge->name]; 

        SIZETYPE sr1 = 0, sr2 = rankvec.bwt_sz;
        simnodes.emplace_back(0, sr1, sr2, FG.size()); // the root SimNode has no base, we simply set it to 0, which does not mean that it has base #(0) 
        for (SIZETYPE w = 0; w <= O.size(); ++w)
        {
            if (w > 0)
                croE[w] = std::max(croE[w - 1] + edge->ue, FG[w - 1].G + edge->ve);
            else
                croE[w] = -inf;
            SCORETYPE tgw = (O.size() > w ? (O.size() - w - 1) * edge->tail->ue + edge->tail->ve : 0);
            FG.emplace_back(w, -inf, std::max(croE[w], Avals[w] + edge->T));
            Ethres[w] = std::max(croE[w], FG[w].G + std::min({SCORETYPE(0), edge->tail->ve - edge->ue, SCORETYPE(tgw - (O.size() - w) * edge->ue)}));

            if (FG[w].G >= edge->head->Avals[0][w])
            {
                if (FG[w].G > edge->head->Avals[0][w])
                {
                    edge->head->Avals[0][w] = FG[w].G;
                    edge->head->AdeltaDot[w].clear();
                    edge->head->AdeltaGlobal[w].clear();
                }
                edge->head->AdeltaGlobal[w].emplace_back(edge, simnodes[0].sr1, 0);
            }
        }

        while (true)
        {
            while (true)
            {
                NUCTYPE c;
                if (simnodes.back().shiftFG < FG.size())
                    c = 1;
                else
                {
                    do
                    {
                        c = simnodes.back().c;
                        if (FG.size() > simnodes.back().shiftFG)
                            FG.erase(FG.begin() + simnodes.back().shiftFG, FG.end());
                        simnodes.pop_back();
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
                    simnodes.emplace_back(c, sr1, sr2, FG.size());
                    break;
                }
                else if (FG.size() > simnodes.back().shiftFG)
                    FG.erase(FG.begin() + simnodes.back().shiftFG, FG.end());
            }
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

                    if (G_val >= edge->head->Avals[0][w])
                    {
                        if (G_val > edge->head->Avals[0][w])
                        {
                            edge->head->Avals[0][w] = G_val;
                            edge->head->AdeltaDot[w].clear();
                            edge->head->AdeltaGlobal[w].clear();
                        }
                        edge->head->AdeltaGlobal[w].emplace_back(edge, simnodes.back().sr1, simnodes.size() - 1);
                    }
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
            Dvals[w - 1] = Avals[w - 1] + edge->T;

        for (SIZETYPE s = 0; s <= ref_sz; ++s)
        {
            if (w > 0)
            {
                SCORETYPE scores[3] = {G0vals[s][w - 1], Dvals[w - 1] + edge->gf[s], -inf};
                if (s > 0)
                    scores[2] = Avals[w - 1] + edge->gfpT[s];
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
                vals[M] = std::max(Evals[s][w - 1] + edge->ue, Gvals[s][w - 1] + edge->ve);
                if (vals[M] == Evals[s][w - 1] + edge->ue)
                    bits[M] += 1;
                if (vals[M] == Gvals[s][w - 1] + edge->ve)
                    bits[M] += 2;
            }
                
            M = &F0vals[s][w] - vals.get();
            bits[M] = 0;
            if (s == 0)
                vals[M] = -inf;
            else
            {
                vals[M] = std::max(F0vals[s - 1][w] + edge->uf, G0vals[s - 1][w] + edge->vf);
                if (vals[M] == F0vals[s - 1][w] + edge->uf)
                    bits[M] += 1;
                if (vals[M] == G0vals[s - 1][w] + edge->vf)
                    bits[M] += 2;
            }

            SCORETYPE scores[3] = {Evals[s][w], F0vals[s][w], -inf};
            if (s > 0 && w > 0)
                scores[2] = Gvals[s - 1][w - 1] + edge->gamma[ref[s - 1]][O[w - 1]];
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

            if (G0vals[s][w] >= edge->head->Avals[0][w])
            {
                if (G0vals[s][w] > edge->head->Avals[0][w])
                {
                    edge->head->Avals[0][w] = G0vals[s][w];
                    edge->head->AdeltaDot[w].clear();
                    edge->head->AdeltaGlobal[w].clear();
                }
                edge->head->AdeltaDot[w].emplace_back(&G0vals[s][w]);
            }
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

            if (circuits[i]->sncs.front().G0 >= circuits[i]->head->Avals[0][w])
            {
                if (circuits[i]->sncs.front().G0 > circuits[i]->head->Avals[0][w])
                {
                    circuits[i]->head->Avals[0][w] = circuits[i]->sncs.front().G0;
                    circuits[i]->head->AdeltaDot[w].clear();
                    circuits[i]->head->AdeltaGlobal[w].clear();
                }
                circuits[i]->head->AdeltaGlobal[w].emplace_back(circuits[i], circuits[i]->sncs.front().sr1, circuits[i]->sncs.front().lambda);
            }

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
                    {
                        if (jumps[i]->G0 >= circuits[i]->head->Avals[0][w])
                        {
                            if (jumps[i]->G0 > circuits[i]->head->Avals[0][w])
                            {
                                circuits[i]->head->Avals[0][w] = jumps[i]->G0;
                                circuits[i]->head->AdeltaDot[w].clear();
                                circuits[i]->head->AdeltaGlobal[w].clear();
                            }
                            circuits[i]->head->AdeltaGlobal[w].emplace_back(circuits[i], jumps[i]->sr1, jumps[i]->lambda);
                        }
                    }
                        
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
            if (file2long.count(edges[swn.n]->name))
                return VALTYPE::LID0;
            if (edges[swn.n]->tail->scc_id==edges[swn.n]->head->scc_id)
            {
                EdgeLocalCircuit *edge = (EdgeLocalCircuit *) edges[swn.n];
                if (M < edge->Evals[0] - vals.get())
                    return VALTYPE::SID;
                if (M < edge->F0vals[0] - vals.get())
                    return VALTYPE::SIE;
                if (M < edge->G0vals[0] - vals.get())
                    return VALTYPE::SIF0;
                if (M < edge->Gvals[0] - vals.get())
                    return VALTYPE::SIG0;
                if (M < edge->D0vals[0]-vals.get())
                    return VALTYPE::SIG;
                if (M < edge->DXvals[0]-vals.get())
                    return VALTYPE::SID0;
                return VALTYPE::SIDX;
            }
            else
            {
                EdgeLocalCross *edge = (EdgeLocalCross *) edges[swn.n];
                if (M < edge->Fvals[0] - vals.get())
                    return VALTYPE::SOE;
                if (M < edge->Gvals[0] - vals.get())
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
                Node &node = nodes[Dot::DOTTYPE2NODEIDX(swn.n)];
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

    struct TrackTree
    {
        DOTTYPE n;
        QUERYSIZE W, lambda = 0;
        std::deque<Dot> dots;
        int64_t idx, pidx;

        struct TrackNode
        {
            int64_t cidxs[5] = {-1, -1, -1, -1, -1};
            QUERYSIZE tau = 0;
        };

        std::deque<TrackNode> tracknodes;

        TrackTree(DOTTYPE n_, QUERYSIZE W_) : n(n_), W(W_)
        {
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

    void GetMinimalGraph()
    {
        std::map<Edge *, TrackTree> long2tracktree;
        std::map<Edge *, std::unique_ptr<SCORETYPE []>> edge2tailAvals;
        for (Edge *edge : edges)
        {
            if (file2long.count(edge->name))
                long2tracktree.try_emplace(edge, edge->n, O.size());
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
                    GlobalTrack(edge2tailAvals, Mval, M, swn, type, globalqueue, localqueue, valuequeue, id, dot_sources, long2tracktree);
                else
                {
                    outputdot(Mval, M, swn, type);
                    switch (type)
                    {
                    case VALTYPE::Q:
                    {
                        SIZETYPE ei = 1;
                        for (Node &node : nodes)
                        {
                            if (!node.is_target)
                                continue;
                            if (ei > std::numeric_limits<BITSTYPE>::max())
                                break;
                            if (bits[M] & ei)
                                arrangesource(valuequeue, &node.Avals[node.scc_sz - 1][O.size()], localqueue, id);
                            ei *= 2;
                        }
                        break;
                    }
                    case VALTYPE::ABAR:
                        break;
                    case VALTYPE::SOE:
                    {
                        EdgeLocalCross *edge = (EdgeLocalCross *)edges[swn.n];
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge->Gvals[swn.s][swn.w - 1]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SOF:
                    {
                        EdgeLocalCross *edge = (EdgeLocalCross *)edges[swn.n];
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-Omax-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge->Gvals[swn.s - 1][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SOG:
                    {
                        EdgeLocalCross *edge = (EdgeLocalCross *)edges[swn.n];
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &(edge->Evals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge->Fvals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & 4)
                            arrangesource(valuequeue, &(edge->tail->Avals[edge->tail->scc_sz - 1][swn.w]), localqueue, id);
                        if (bits[M] & 8)
                            arrangesource(valuequeue, &vals[M-Omax-2], localqueue, id);
                        break;
                    }
                    case VALTYPE::SID:
                        arrangesource(valuequeue, &(edges[swn.n]->tail->Avals[edges[swn.n]->tail->scc_sz - 1][swn.w]), localqueue, id);
                        break;
                    case VALTYPE::SIE:
                    {
                        EdgeLocalCircuit *edge = (EdgeLocalCircuit *)edges[swn.n];
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge->Gvals[swn.s][swn.w - 1]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIF0:
                    {
                        EdgeLocalCircuit *edge = (EdgeLocalCircuit *)edges[swn.n];
                        if (bits[M] & 1)
                            arrangesource(valuequeue, &vals[M-Omax-1], localqueue, id);
                        if (bits[M] & 2)
                            arrangesource(valuequeue, &(edge->G0vals[swn.s - 1][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIG0:
                    {
                        EdgeLocalCircuit *edge = (EdgeLocalCircuit *)edges[swn.n];
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &(edge->Evals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(edge->F0vals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(4))
                            arrangesource(valuequeue, &(edge->Gvals[swn.s - 1][swn.w - 1]), localqueue, id);
                        break;
                    }
                    case VALTYPE::SIG:
                    {
                        EdgeLocalCircuit *edge = (EdgeLocalCircuit *)edges[swn.n];
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &(edge->G0vals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(edge->Dvals[swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(4))
                            arrangesource(valuequeue, &(edge->tail->Avals[edge->tail->scc_sz - 1][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::LID0: case VALTYPE::SID0:
                        arrangesource(valuequeue, &(edges[swn.n]->tail->Avals[swn.s][swn.w]), localqueue, id);
                        break;
                    case VALTYPE::SIDX:
                    {
                        EdgeLocalCircuit *edge = (EdgeLocalCircuit *) edges[swn.n];
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &(edge->tail->Avals[swn.s][swn.w]), localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(edge->D0vals[swn.s][swn.w]), localqueue, id);
                        break;
                    }
                    case VALTYPE::NA:
                    {
                        Node &node = nodes[Dot::DOTTYPE2NODEIDX(swn.n)];
                        if (bits[M] & int8_t(1))
                            arrangesource(valuequeue, &node.Avals[0][swn.w], localqueue, id);
                        SIZETYPE ei = 1;
                        for (Edge *edge : edges)
                        {
                            if (edge->head != &node || edge->tail->scc_id != node.scc_id)
                                continue;
                            ei *= 2;
                            if (ei > std::numeric_limits<BITSTYPE>::max())
                                break;
                            if (bits[M] & int8_t(ei))
                            {
                                SCORETYPE *source;
                                if (file2short.count(edge->name))
                                {
                                    EdgeLocalCircuit *circuit = (EdgeLocalCircuit *)edge;
                                    if (Mval == circuit->D0vals[swn.s - 1][swn.w] + circuit->gfm) // this works because D0/DX[l-1] are only tracked by A[l]
                                        source = &(circuit->D0vals[swn.s - 1][swn.w]);
                                    else
                                        source = &(circuit->DXvals[swn.s - 1][swn.w]);
                                }
                                else
                                {
                                    EdgeGlobalCircuit *circuit = (EdgeGlobalCircuit *)edge;
                                    source = &(circuit->D0vals[swn.s - 1][swn.w]);
                                }
                                arrangesource(valuequeue, source, localqueue, id);
                            }
                        }
                        break;
                    }
                    case VALTYPE::NB:
                    {
                        DOTTYPE n = Dot::DOTTYPE2NODEIDX(swn.n);
                        if (bits[M] & uint8_t(1))
                            arrangesource(valuequeue, &vals[M-1], localqueue, id);
                        if (bits[M] & uint8_t(2))
                            arrangesource(valuequeue, &(nodes[n].Avals[nodes[n].scc_sz - 1][swn.w - 1]), localqueue, id);
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
                        Node *tail = edges[M->n]->tail;
                        arrangesource(valuequeue, &tail->Avals[tail->scc_sz - 1][M->w], localqueue, id);
                    }
                }
            }
        }
        if (id > max_id)
            max_id = id;

        for (Node &node : nodes)
            for (SIZETYPE w = 0; w <= O.size(); ++w)
            {
                node.AdeltaDot[w].clear();
                node.AdeltaGlobal[w].clear();
            }
    }

    void GlobalTrack(std::map<Edge *, std::unique_ptr<SCORETYPE []>> &edge2tailAvals, SCORETYPE Mval, SIZETYPE M, SWN &swn, Align::VALTYPE type, std::queue<Dot *> &globalqueue, std::queue<SIZETYPE> &localqueue, std::queue<SCORETYPE> &valuequeue, IDTYPE &id, std::deque<Dot *> &dot_sources, std::map<Edge *, TrackTree> &long2tracktree)
    {
        Node &node = nodes[Dot::DOTTYPE2NODEIDX(swn.n)];
        outputdot(Mval, M, swn, type);
        for (SCORETYPE *source : node.AdeltaDot[swn.w])
            arrangesource(valuequeue, source, localqueue, id);

        for (Node::GlobalSuffix &globalsuffix : node.AdeltaGlobal[swn.w])
        {
            Edge *edge = globalsuffix.edge;
            SCORETYPE *Avals = edge2tailAvals[edge].get();
            TrackTree &tracktree = long2tracktree.at(edge);
            std::deque<TrackTree::TrackNode> &tracknodes = tracktree.tracknodes;
            RankVec &rankvec = file2rankvec[edge->name];
            std::ifstream &SAfin = file2SA[edge->name];
            SIZETYPE start;
            SAfin.seekg(globalsuffix.start * sizeof(SIZETYPE));
            SAfin.read((char*)&start, sizeof(SIZETYPE));
            std::ifstream &REVREFfin = file2long[edge->name];
            NUCTYPE revref[globalsuffix.lambda];
            REVREFfin.seekg(start);
            REVREFfin.read((char*)revref, globalsuffix.lambda);

            tracktree.pidx = -1;
            tracktree.idx = 0;
            tracktree.lambda = 0;
            if (edge->head->scc_id != edge->tail->scc_id)
            {
                for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w == 0)
                    {
                        tracktree.dotE(w).val = -inf;
                        tracktree.dotE(w).s_sz = 0;
                    }
                    else
                        source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge->ue, tracktree.dotG(w - 1).val + edge->ve}, dot_sources);
                    tracktree.dotF(w).val = -inf;
                    tracktree.dotF(w).s_sz = 0;
                    source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w), NULL}, {tracktree.dotE(w).val, tracktree.dotF(w).val, Avals[w] + edge->T}, dot_sources);
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
                            source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge->ue, tracktree.dotG(w - 1).val + edge->ve}, dot_sources);
                        source_max(tracktree.dotF(w), {&tracktree.dotPF(w), &tracktree.dotPG(w)}, {tracktree.dotPF(w).val + edge->uf, tracktree.dotPG(w).val + edge->vf}, dot_sources);
                        if (w == 0)
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w)}, {tracktree.dotE(w).val, tracktree.dotF(w).val}, dot_sources);
                        else
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w), &tracktree.dotPG(w - 1)}, {tracktree.dotE(w).val, tracktree.dotF(w).val, tracktree.dotPG(w - 1).val + edge->gamma[c][O[w - 1]]}, dot_sources);
                    }
                    tracknodes[tracktree.idx].tau = std::max(tracknodes[tracktree.idx].tau, swn.w + 1);
                }
            }
            else
            {
                for (SIZETYPE w = tracknodes[tracktree.idx].tau; w <= swn.w; ++w)
                {
                    if (w > 0)
                        source_max(tracktree.dotG(w - 1), {&tracktree.dotE(w - 1), NULL}, {tracktree.dotE(w - 1).val, Avals[w - 1] + edge->T}, dot_sources);
                    if (w == 0)
                    {
                        tracktree.dotE(w).val = -inf;
                        tracktree.dotE(w).s_sz = 0;
                    }
                    else
                        source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge->ue, tracktree.dotG(w - 1).val + edge->ve}, dot_sources);
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
                            source_max(tracktree.dotE(w), {&tracktree.dotE(w - 1), &tracktree.dotG(w - 1)}, {tracktree.dotE(w - 1).val + edge->ue, tracktree.dotG(w - 1).val + edge->ve}, dot_sources);
                        if (tracktree.pidx == 0)
                            source_max(tracktree.dotF(w), {&tracktree.dotPF(w), &tracktree.dotPE(w)}, {tracktree.dotPF(w).val + edge->uf, tracktree.dotPE(w).val + edge->vf}, dot_sources);
                        else
                            source_max(tracktree.dotF(w), {&tracktree.dotPF(w), &tracktree.dotPG(w)}, {tracktree.dotPF(w).val + edge->uf, tracktree.dotPG(w).val + edge->vf}, dot_sources);
                        if (w == 0)
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w)}, {tracktree.dotE(w).val, tracktree.dotF(w).val}, dot_sources);
                        else
                            source_max(tracktree.dotG(w), {&tracktree.dotE(w), &tracktree.dotF(w), &tracktree.dotPG(w - 1)}, {tracktree.dotE(w).val, tracktree.dotF(w).val, tracktree.dotPG(w - 1).val + edge->gamma[c][O[w - 1]]}, dot_sources);
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
