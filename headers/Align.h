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
struct Truther
{
    std::istringstream & GetDivideLine(std::ifstream & fin, std::string & line, std::istringstream & iss)
    {
        std::getline(fin,line);
        iss.clear();
        iss.str(line);
        return iss;
    }

    void PathToTruth(std::string file, std::string out_file, int max_glo)
    {
        std::ifstream fin(file);
        std::ofstream fout(out_file);
        std::string line;
        std::istringstream iss;
        T* glo_poses=new T[max_glo];
        int eid,s,w,st,pres,now_glo;
        double val;
        getline(fin,line);
        while(fin.good())
        {
            fout << line << '\n';
            do
            {
                GetDivideLine(fin, line, iss) >> eid >> s >> w >> val;
            }while(fin.good() && eid==-1);
            while(!iss.fail())
            {
                fout << eid;
                st=s;
                if(iss >> now_glo)
                    for(int i=0; i<now_glo; ++i)
                        iss >> glo_poses[i];
                else
                    now_glo=1, glo_poses[0]=0;
                do
                {
                    pres=s;
                    GetDivideLine(fin, line, iss) >> eid >> s;
                }while(s>=0);
                for(int i=0; i<now_glo; ++i)
                    fout << '\t' << st+glo_poses[i] << '\t' << pres+glo_poses[i];
                fout << '\n';
                do
                {
                    GetDivideLine(fin, line, iss) >> eid >> s >> w >> val;
                }while(fin.good() && (eid==-1 || s<0));
            }
            iss.clear();
        }
        fin.close();
        delete[] glo_poses;
    }
};

template <typename T>
struct Tracker
{
    std::tuple<int,int*,Dot<T>*,T*,Dot<T>**> ReadMinimalGraph(std::ifstream & fin)
    {
        int nt;
        fin.read((char*)&nt,sizeof(int));
        if(!fin.good())
            return std::tuple<int,int*,Dot<T>*,T*,Dot<T>**>();
        int* tids=new int[nt];
        fin.read((char*)tids,nt*sizeof(int));
        size_t nd;
        fin.read((char*)&nd,sizeof(size_t));
        Dot<T>* VV=new Dot<T>[nd];
        fin.read((char*)VV,nd*sizeof(Dot<T>));
        int robs_sz;
        fin.read((char*)&robs_sz,sizeof(int));
        T** robs=new T*[robs_sz];
        fin.read((char*)robs,robs_sz*sizeof(T*));
        T* ranges=new T[3*robs_sz];
        fin.read((char*)ranges,3*robs_sz*sizeof(T));
        std::map<T*,T*> robs2ranges;
        for(int i=0; i<robs_sz; ++i)
            robs2ranges.emplace(robs[i],ranges+3*i);
        int id_sz;
        fin.read((char*)&id_sz,sizeof(int));
        int* ids=new int[id_sz];
        fin.read((char*)ids,id_sz*sizeof(int));
        Dot<T>** sources=new Dot<T>*[id_sz];
        int idd=0;
        Dot<T>** sources_now=sources;
        for(size_t i=0; i<nd; ++i)
        {
            Dot<T> & dot=VV[i];
            if(dot.rob)
                dot.rob=robs2ranges[dot.rob];
            dot.sources=sources_now;
            for(int j=0; j<dot.s_sz; ++j,++idd)
                sources_now[j]=&VV[ids[idd]];
            sources_now+=dot.s_sz;
        }
        delete[] ids;
        return std::make_tuple(nt,tids,VV,ranges,sources);
    }

    void ExplorePath(Dot<T>* target, int read_num, int target_num, int max_num, int max_glo, Graph<T> & graph, std::ofstream & fout)
    {
        std::vector<Dot<T>*> PP;
        std::stack<Dot<T>*> TT;
        TT.push(target);
        int accum=0;
        while(!TT.empty())
        {
            PP.push_back(TT.top());
            TT.pop();
            if(PP.back()->s_sz>0)
                for(int i=PP.back()->s_sz-1; i>=0; --i)
                    TT.push(PP.back()->sources[i]);
            else
            {
                fout << read_num << '\t' << target_num << '\t' << accum << '\n';
                for(auto it=PP.rbegin(); it!=PP.rend(); ++it)
                {
                    fout << (*it)->eid << '\t' << (*it)->s << '\t' << (*it)->w << '\t' << (*it)->val;
                    if((*it)->rob && (*(it-1))->eid==-1)
                    {
                        int now_glo=std::min(int((*it)->rob[1]-(*it)->rob[0]+1),max_glo);
                        fout << '\t' << now_glo;
                        BroWheel<T> & Bro=graph.eid2globals[(*it)->eid]->Bro;
                        for(int i=(*it)->rob[0]; i<(*it)->rob[0]+now_glo; ++i)
                        {
                            T spos=Bro.size()-1-Bro.SimSuffix(i)-(*it)->rob[2];
                            fout << '\t' << spos;
                        }
                    }
                    fout << '\n';
                }
                ++accum;
                if(accum==max_num)
                    break;
                Dot<T>* ptr;
                do
                {
                    ptr=PP.back();
                    PP.pop_back();
                }
                while(!PP.empty() && ptr==PP.back()->sources[PP.back()->s_sz-1]);
            }
        }
    }

    void ReadAndExploreAll(std::string file, std::string out_file, int max_num, int max_glo, Graph<T> & graph)
    {
        std::ifstream fin(file,std::ifstream::binary);
        std::ofstream fout(out_file);
        fout << std::setprecision(10);
        int nt, *tids;
        Dot<T> *VV, **sources;
        T* ranges;
        std::tie(nt,tids,VV,ranges,sources)=ReadMinimalGraph(fin);
        for(int read_num=0; fin.good(); ++read_num)
        {
            for(int target_num=0; target_num<nt; ++target_num)
            {
                ExplorePath(&VV[tids[target_num]], read_num, target_num, max_num, max_glo, graph, fout);
            }
            delete[] tids;
            delete[] VV;
            delete[] ranges;
            delete[] sources;
            std::tie(nt,tids,VV,ranges,sources)=ReadMinimalGraph(fin);
        }
        fin.close();
        fout.close();
    }
};

struct Do
{
    int w;
    double val;
};

template <typename T>
struct TrieState
{
    Do* E;
    Do* F;
    Do* G;
    int Efn;
    int Ffn;
    int Gfn;
    char* now;
    size_t remain;
    size_t chunk_id;
};

template <typename T>
struct SimNode
{
    bool visit;
    T letter;
    T sr1;
    T sr2;
};

template <typename U>
struct Align : Memory, GraphCopy<U>
{
    Graph<U> & graph;
    std::vector<Dot<U>*> adrs;
    std::vector<double> vals;
    std::vector<Dot<U>*> VV;
    std::ofstream MGout;
    std::vector<TrieState<U>> triestates;
    std::stack<SimNode<U>> simnodes;
    
    Align(Graph<U> & graph_v) : graph(graph_v)
    {
        this->GetCopy(graph_v);
        Initial(1024*1024*1024);
    }
    
    void CallGetMinimalGraph(int W,std::string & O)
    {
        MGout.write((char*)&(graph.t_sz),sizeof(int));
        VV.clear();
        for(int i=0; i<graph.t_sz; ++i)
        {
            Node<U> & node=graph.nodes[graph.targets[i]];
            Dot<U>* target=&(this->nodecopys[graph.targets[i]].tildeA[node.tAsz-1][W]);
            if(!target->visit)
                GetMinimalGraph(target,O);
            MGout.write((char*)&(target->id),sizeof(int));
        }
        size_t nd=VV.size();
        MGout.write((char*)&(nd),sizeof(size_t));
        std::set<U*> robs;
        int id_sz=0;
        for(Dot<U>* dot : VV)
        {
            MGout.write((char*)dot,sizeof(Dot<U>));
            if(dot->rob)
                robs.insert(dot->rob);
            id_sz+=dot->s_sz;
        }
        int robs_sz=robs.size();
        MGout.write((char*)&robs_sz,sizeof(int));
        for(U* rob : robs)
            MGout.write((char*)&rob,sizeof(U*));
        for(U* rob : robs)
            MGout.write((char*)rob,3*sizeof(U));
        MGout.write((char*)&id_sz,sizeof(int));
        for(Dot<U>* dot : VV)
        {
            for(int i=0; i<dot->s_sz; ++i)
                MGout.write((char*)&(dot->sources[i]->id),sizeof(int));
        }
    }
    
    void GetMinimalGraph(Dot<U>* target, std::string & O)
    {
        target->visit=true;
        GlobalTrack(target,O);
        for(int i=0; i<target->s_sz; ++i)
            if(!target->sources[i]->visit)
                GetMinimalGraph(target->sources[i],O);
        VV.push_back(target);
        target->id=VV.size()-1;
    }
    
    void GlobalTrack(Dot<U>* target, std::string & O)
    {
        if(target->s!=-1 || target->eid==-1)
            return;
        EdgeGlobal<U>* global=graph.eid2globals[target->eid];
        if(!global)
            return;
        auto & box=this->globalcopys[global-graph.globals].boxes[target->w];
        if(box.empty())
            return;
        Node<U> & NT=graph.nodes[global->tail];
        NodeCopy<U> & NTCP=this->nodecopys[global->tail];
        for(size_t i=0; i<box.size(); i+=3)
        {
            std::string seq;
            seq.resize(box[i+2]);
            U os=global->Bro.SimSuffix(box[i]);
            adrs.clear();
            for(U j=os,k=box[i+2]-1; j<os+box[i+2]; ++j,--k)
                seq[k]=global->G(j);
            U* rob=heap_alloc<U>(3);
            rob[0]=box[i];
            rob[1]=box[i+1];
            rob[2]=box[i+2];
            Dot<U>** E=heap_alloc<Dot<U>*>(rob[2]+1);
            Dot<U>** F=heap_alloc<Dot<U>*>(rob[2]+1);
            Dot<U>** G=heap_alloc<Dot<U>*>(rob[2]+1);
            alloc_initial(E, target->eid, rob[2], target->w, rob);
            alloc_initial(F, target->eid, rob[2], target->w, rob);
            alloc_initial(G, target->eid, rob[2], target->w, rob);
            if(graph.cross[target->eid])
            {
                CrossInitial(E,F,G,NTCP.tildeA[NT.tAsz-1],target->w,global->ve,global->ue,global->T);
                for(int s=1; s<=rob[2]; ++s)
                    CrossBody(seq,O,s,target->w,E,F,G,global->ve,global->ue,NULL,NULL,*global);
            }
            else
            {
                Dot<U>** H2=heap_alloc<Dot<U>*>(1);
                alloc_initial(H2, target->eid, 0, target->w, rob);
                H2[0][target->w].val=-inf;
                for(size_t s=0; s<=seq.size(); ++s)
                    E[s][0].val=F[s][0].val=G[s][0].val=-inf;
                for(int wp=1; wp<=target->w; ++wp)
                {
                    source_max(H2[0][wp-1],{&NTCP.tildeA[NT.tAsz-1][wp-1]},{NTCP.tildeA[NT.tAsz-1][wp-1].val+global->T});
                    CircleInitial(wp,E,F,G,H2[0][wp-1],global->ve,global->ue);
                    CircleBody(seq,O,wp,1,E,F,G,H2[0][wp-1],H2[0][target->w],global->ve,global->ue,*global);
                    for(size_t s=2; s<=seq.size(); ++s)
                        CircleBody(seq,O,wp,s,E,F,G,H2[0][target->w],H2[0][target->w],global->ve,global->ue,*global);
                }
            }
            adrs.push_back(&G[rob[2]][target->w]);
        }
        Dot<U>** sources=target->sources;
        int s_sz=target->s_sz;
        target->s_sz+=adrs.size();
        target->sources=heap_alloc<Dot<U>*>(target->s_sz);
        int j=0;
        for(int i=0; i<s_sz; ++i,++j)
            target->sources[j]=sources[i];
        for(size_t i=0; i<adrs.size(); ++i,++j)
            target->sources[j]=adrs[i];
        box.clear();
    }
    
    void Mix(std::string & O)
    {
        clear();
        for(int i=0; i<graph.n_sz; ++i)
        {
            Node<U> & node=graph.nodes[i];
            if(node.tAsz>0)
            {
                for(int j=0; j<node.tAsz; ++j)
                    this->nodecopys[i].tildeA[j]=alloc_initial(-1, O.size());
                this->nodecopys[i].cid_now=0;
                size_t ge_sz=node.ge.size();
                if(ge_sz<O.size()+1)
                {
                    node.ge.resize(O.size()+1);
                    for(size_t j=ge_sz; j<=O.size(); ++j)
                        node.ge[j]=node.ge[j-1]+node.tue;
                }
            }
        }
        for(int i=0; i<graph.r_sz; ++i)
        {
            Dot<U>* barA=this->nodecopys[graph.roots[i]].barA=alloc_initial(-1, O.size());
            for(size_t w=0; w<=O.size(); ++w)
            {
                barA[w].val=graph.nodes[graph.roots[i]].ge[w];
                barA[w].s_sz=0;
            }
        }
        for(int i=0; i<graph.n_sz; ++i)
            update_cross(this->nodecopys[i], O.size());
        for(int i=0,l=0,g=0; i<graph.scc_n; ++i)
        {
            for(int j=0; j<graph.local_C[i]; ++j,++l)
            {
                EdgeLocal<U> & edge=graph.locals[l];
                EdgeLocalCopy<U> & edgecopy=this->localcopys[l];
                EdgeAllocSpace(edge,edgecopy,O.size(),graph.tAsz[i]);
                alloc_initial(edgecopy.E,edge.id,edge.seq.size(),O.size(),NULL);
                alloc_initial(edgecopy.F,edge.id,edge.seq.size(),O.size(),NULL);
                alloc_initial(edgecopy.G,edge.id,edge.seq.size(),O.size(),NULL);
            }
            for(int j=0; j<graph.global_C[i]; ++j,++g)
            {
                EdgeGlobal<U> & edge=graph.globals[g];
                EdgeGlobalCopy<U> & edgecopy=this->globalcopys[g];
                EdgeAllocSpace(edge,edgecopy,O.size(),graph.tAsz[i]);
                if(edgecopy.boxes.size()<O.size()+1)
                    edgecopy.boxes.resize(O.size()+1);
            }
        }
        for(int i=0,l=0,g=0; i<graph.scc_n; ++i)
        {
            if(graph.local_C[i]==1 && graph.global_C[i]==0 && graph.locals[l].tail!=graph.locals[l].head)
                CrossIteration(l++, O);
            else
                if(graph.local_C[i]==0 && graph.global_C[i]==1 && graph.globals[g].tail!=graph.globals[g].head)
                    CrossIterationGlobal(g++, O);
                else
                {
                    CircleIteration(l,g,i,O);
                    l+=graph.local_C[i];
                    g+=graph.global_C[i];
                }
        }
    }
    
    template <typename TP, typename TPC>
    void EdgeAllocSpace(TP & edge, TPC & edgecopy, size_t W, int tAsz)
    {
        if(edge.vf.size()<W+1)
        {
            edge.vf.resize(W+1,edge.vf.back());
            edge.uf.resize(W+1,edge.uf.back());
        }
        for(int a=0; a<tAsz; ++a)
            edgecopy.A[a]=alloc_initial(edge.id, W, a==0?-1:-2);
        edgecopy.B=alloc_initial(edge.id, W,-3);
    }
    
    void CrossIteration(int l, std::string & O)
    {
        EdgeLocal<U> & edge=graph.locals[l];
        EdgeLocalCopy<U> & edgecopy=this->localcopys[l];
        Dot<U>* tildeANT=this->nodecopys[edge.tail].tildeA[graph.nodes[edge.tail].tAsz-1];
        if(tildeANT[0].val==-inf)
            return;
        
        CrossInitial(edgecopy.E,edgecopy.F,edgecopy.G,tildeANT,O.size(),edge.ve[0],edge.ue[0],edge.T);
        for(size_t s=1; s<=edge.seq.size(); ++s)
        {
            CrossBody(edge.seq,O,s,O.size(),edgecopy.E,edgecopy.F,edgecopy.G,edge.ve[s],edge.ue[s],tildeANT,edge.gfp,edge);
        }
        edgecopy.B[0].val=-inf;
        adrs.clear();
        vals.clear();
        for(size_t s=0; s<=edge.seq.size(); ++s)
        {
            adrs.push_back(&edgecopy.G[s][0]);
            vals.push_back(edgecopy.G[s][0].val+edge.gfm[s]);
        }
        source_max(edgecopy.A[0][0],adrs.begin(),vals.begin(),vals.end());
        for(size_t w=1; w<=O.size(); ++w)
        {
            source_max(edgecopy.B[w],{&edgecopy.B[w-1],&edgecopy.A[0][w-1]},{edgecopy.B[w-1].val+graph.nodes[edge.head].tue,edgecopy.A[0][w-1].val+graph.nodes[edge.head].tve});
            adrs.clear();
            vals.clear();
            adrs.push_back(&edgecopy.B[w]);
            vals.push_back(edgecopy.B[w].val);
            for(size_t s=0; s<=edge.seq.size(); ++s)
            {
                adrs.push_back(&edgecopy.G[s][w]);
                vals.push_back(edgecopy.G[s][w].val+edge.gfm[s]);
            }
            source_max(edgecopy.A[0][w],adrs.begin(),vals.begin(),vals.end());
        }
        ++this->nodecopys[edge.head].cid_now;
        update_cross(this->nodecopys[edge.head], O.size());
    }
    
    void CrossInitial(Dot<U>** E, Dot<U>** F, Dot<U>** G, Dot<U>* tildeANT, int W, double ve, double ue, double T)
    {
        for(int j=0; j<=W; ++j)
            F[0][j].val=-inf;
        E[0][0].val=-inf;
        source_max(G[0][0],{&tildeANT[0]},{tildeANT[0].val+T});
        for(int w=1; w<=W; ++w)
        {
            source_max(E[0][w],{&E[0][w-1],&G[0][w-1]},{E[0][w-1].val+ue,G[0][w-1].val+ve});
            source_max(G[0][w],{&E[0][w],&tildeANT[w]},{E[0][w].val,tildeANT[w].val+T});
        }
    }
    
    template <typename Y>
    void CrossBody(std::string & seq, std::string & O, int s, int W, Dot<U>** E, Dot<U>** F, Dot<U>** G, double ve, double ue, Dot<U>* tildeANT, double* gfp, Y & edge)
    {
        for(int w=0; w<=W; ++w)
            source_max(F[s][w],{&F[s-1][w],&G[s-1][w]},{F[s-1][w].val+edge.uf[w],G[s-1][w].val+edge.vf[w]});
        E[s][0].val=-inf;
        if(tildeANT)
            source_max(G[s][0],{&F[s][0],&tildeANT[0]},{F[s][0].val,tildeANT[0].val+gfp[s]+edge.T});
        else
            source_max(G[s][0],{&F[s][0]},{F[s][0].val});
        for(int w=1; w<=W; ++w)
        {
            source_max(E[s][w],{&E[s][w-1],&G[s][w-1]},{E[s][w-1].val+ue,G[s][w-1].val+ve});
            if(tildeANT)
                source_max(G[s][w],{&F[s][w],&E[s][w],&G[s-1][w-1],&tildeANT[w]},{F[s][w].val,E[s][w].val,G[s-1][w-1].val+edge.gamma[int(seq[s-1])][int(O[w-1])],tildeANT[w].val+gfp[s]+edge.T});
            else
                source_max(G[s][w],{&F[s][w],&E[s][w],&G[s-1][w-1]},{F[s][w].val,E[s][w].val,G[s-1][w-1].val+edge.gamma[int(seq[s-1])][int(O[w-1])]});
        }
    }
    
    void CrossIterationGlobal(int g, std::string & O)
    {
        EdgeGlobal<U> & edge=graph.globals[g];
        EdgeGlobalCopy<U> & edgecopy=this->globalcopys[g];
        Dot<U>* tildeANT=this->nodecopys[edge.tail].tildeA[graph.nodes[edge.tail].tAsz-1];
        if(tildeANT[0].val==-inf)
            return;
        Node<U> & NH=graph.nodes[edge.head];
        double tve=NH.tve, tue=NH.tue;
        std::vector<double> & ge=NH.ge;
        BroWheel<U> & Bro=edge.Bro;
        double ve=edge.ve, ue=edge.ue, T=edge.T;
        std::vector<double> &vf=edge.vf, &uf=edge.uf;
        std::vector<std::vector<U>> &boxes=edgecopy.boxes;
        auto & gamma=edge.gamma;
        double gamma_max=graph.gamma_max;
        
        Dot<U> *A=edgecopy.A[0], *B=edgecopy.B;
        for(size_t w=0; w<=O.size(); ++w)
        {
            A[w].val=-inf;
            boxes[w].resize(3);
            boxes[w][0]=boxes[w][2]=0;
            boxes[w][1]=Bro.size()-1;
        }
        simnodes.emplace();
        simnodes.top().sr1=0;
        simnodes.top().sr2=Bro.size()-1;
        
        TrieState<U> & ts=PushTrieState(O.size(),0,O.size()+1);
        double opt=-inf;
        ts.G[0].w=0;
        ts.G[0].val=tildeANT[0].val+T;
        A[0].val=ts.G[0].val;
        opt=ts.G[0].val+ge[O.size()];
        for(size_t g=1,e=0; g<=O.size(); ++g,++e)
        {
            ts.G[g].w=ts.E[e].w=g;
            ts.E[e].val=ts.G[e].val+ve;
            if(e>0)
                ts.E[e].val=std::max(ts.E[e].val,ts.E[e-1].val+ue);
            ts.G[g].val=std::max(ts.E[e].val,tildeANT[g].val+T);
            A[g].val=ts.G[g].val;
            opt=std::max(opt,ts.G[g].val+ge[O.size()-g]);
        }
        simnodes.top().visit=true;
        PushPre(edge);
        
        while(!simnodes.empty())
        {
            if(simnodes.top().visit)
            {
                TrieState<U> & ts=triestates.back();
                this->now=ts.now;
                this->remain=ts.remain;
                this->chunk_id=ts.chunk_id;
                triestates.pop_back();
                simnodes.pop();
            }
            else
            {
                TrieState<U> & tsp=triestates.back();
                int Efn, Ffn=tsp.Gfn, Gfn=0;
                if(Ffn>0)
                    Gfn=O.size()-tsp.G[0].w+1;
                Efn=std::max(0,Gfn-1);
                TrieState<U> & ts=PushTrieState(Efn,Ffn,Gfn);
                for(int g=0,f=0; g<Gfn; ++g)
                {
                    int w=ts.F[g].w=tsp.G[g].w;
                    ts.F[g].val=tsp.G[g].val+vf[w];
                    if(f<tsp.Ffn && w==tsp.F[f].w)
                    {
                        ts.F[g].val=std::max(ts.F[g].val,tsp.F[f].val+uf[w]);
                        ++f;
                    }
                }
                U letter=simnodes.top().letter;
                if(Gfn>0)
                {
                    ts.G[0].w=O.size()-Gfn+1;
                    ts.G[0].val=ts.F[0].val;
                    for(int g=1,f=1,e=0; g<Gfn; ++g,++e)
                    {
                        int w=ts.E[e].w=ts.G[g].w=ts.G[e].w+1;
                        ts.E[e].val=ts.G[e].val+ve;
                        if(e>0)
                            ts.E[e].val=std::max(ts.E[e].val,ts.E[e-1].val+ue);
                        ts.G[g].val=ts.E[e].val;
                        if(tsp.G[f-1].w==w-1)
                            ts.G[g].val=std::max(ts.G[g].val,tsp.G[f-1].val+gamma[letter][int(O[w-1])]);
                        if(f<Ffn && ts.F[f].w==w)
                        {
                            ts.G[g].val=std::max(ts.G[g].val,ts.F[f].val);
                            ++f;
                        }
                    }
                }
                SimNode<U> & sn=simnodes.top();
                for(int g=0; g<Gfn; ++g)
                {
                    int w=ts.G[g].w;
                    if(ts.G[g].val>=A[w].val)
                    {
                        if(ts.G[g].val>A[w].val)
                        {
                            A[w].val=ts.G[g].val;
                            boxes[w].clear();
                            opt=std::max(opt,A[w].val+ge[O.size()-w]);
                        }
                        boxes[w].push_back(sn.sr1);
                        boxes[w].push_back(sn.sr2);
                        boxes[w].push_back(triestates.size()-1);
                    }
                }
                if(Gfn>0)
                {
                    int w=ts.G[0].w;
                    if(ts.G[0].val<tildeANT[w].val+T || ts.G[0].val+(O.size()-w)*gamma_max<opt)
                        ts.G[0].val=ts.F[0].val=-inf;
                    for(int g=1,e=0,f=1;g<Gfn;++g,++e)
                    {
                        int w=ts.G[g].w;
                        if(ts.G[g].val+(O.size()-w)*gamma_max<opt)
                        {
                            ts.G[g].val=ts.E[e].val=-inf;
                            if(f<Ffn && w==ts.F[f].w)
                                ts.F[f].val=-inf;
                        }
                        else
                        {
                            if(ts.G[g].val<tildeANT[w].val+T)
                            {
                                ts.G[g].val=-inf;
                                if(f<Ffn && w==ts.F[f].w)
                                    ts.F[f].val=-inf;
                            }
                            if(ts.E[w].val+ue<tildeANT[w].val+T+ve)
                                ts.E[w].val=-inf;
                        }
                        if(f<Ffn && w==ts.F[f].w)
                            ++f;
                    }    
                }
                {
                    int* Xfns[3]={&ts.Efn,&ts.Ffn,&ts.Gfn};
                    int i=0;
                    for(Do* X : {ts.E,ts.F,ts.G})
                    {
                        int XfnNew=0;
                        for(int x=0; x<*Xfns[i]; ++x)
                            if(X[x].val!=-inf)
                            {
                                if(x!=XfnNew)
                                    X[XfnNew]=X[x];
                                ++XfnNew;
                            }
                        *Xfns[i]=XfnNew;
                        ++i;
                    }
                }
                simnodes.top().visit=true;
                if(ts.Gfn>0)
                    PushPre(edge);
            }
        }
        B[0].val=-inf;
        A[0].s_sz=0;
        for(size_t w=1; w<=O.size(); ++w)
        {
            source_max(B[w],{&B[w-1],&A[w-1]},{B[w-1].val+tue,A[w-1].val+tve});
            if(B[w].val>=A[w].val)
            {
                if(B[w].val>A[w].val)
                    boxes[w].clear();
                source_max(A[w],{&B[w]},{B[w].val});
            }
            else
                A[w].s_sz=0;
        }
        ++this->nodecopys[edge.head].cid_now;
        update_cross(this->nodecopys[edge.head], O.size());
    }
    
    TrieState<U> & PushTrieState(int Efn, int Ffn, int Gfn)
    {
        triestates.emplace_back();
        TrieState<U> & ts=triestates.back();
        ts.now=this->now;
        ts.remain=this->remain;
        ts.chunk_id=this->chunk_id;
        ts.Efn=Efn, ts.Ffn=Ffn, ts.Gfn=Gfn;
        ts.E=heap_alloc<Do>(ts.Efn);
        ts.F=heap_alloc<Do>(ts.Ffn);
        ts.G=heap_alloc<Do>(ts.Gfn);
        return ts;
    }
    
    void PushPre(EdgeGlobal<U> & edge)
    {
        U sr1=simnodes.top().sr1, sr2=simnodes.top().sr2;
        for(U letter=2; letter<7; ++letter)
        {
            U sr1cp=sr1, sr2cp=sr2;
            edge.Bro.PreRange(sr1cp, sr2cp, letter);
            if(sr1cp>sr2cp)
                continue;
            simnodes.emplace();
            simnodes.top().letter=letter;
            simnodes.top().sr1=sr1cp;
            simnodes.top().sr2=sr2cp;
            simnodes.top().visit=false;
        }
    }
    
    void CircleIteration(int l, int g, int scc, std::string & O)
    {
        int tAsz=graph.tAsz[scc];
        {
            int i;
            for(i=0; i<graph.n_sz; ++i)
                if(graph.nodes[i].scc==scc && this->nodecopys[i].tildeA[0][0].val!=-inf)
                    break;
            if(i==graph.n_sz)
            {
                for(int j=0; j<graph.n_sz; ++j)
                    for(size_t w=0; w<=O.size(); ++w)
                        this->nodecopys[j].tildeA[tAsz-1][w].val=-inf;
                return;
            }
        }
        for(int ll=l; ll<l+graph.local_C[scc]; ++ll)
        {
            EdgeLocal<U> & edge=graph.locals[ll];
            EdgeLocalCopy<U> & edgecopy=this->localcopys[ll];
            for(int a=0; a<tAsz-1; ++a)
            {
                edgecopy.H1a[a]=alloc_initial(edge.id,O.size(),0);
                edgecopy.H1b[a]=alloc_initial(edge.id,O.size(),edge.seq.size());
            }
            alloc_initial(edgecopy.H2,edge.id,edge.seq.size(),O.size(),NULL);
            alloc_initial(edgecopy.C,edge.id,edge.seq.size(),O.size(),NULL);
        }
        for(int gg=g; gg<g+graph.global_C[scc]; ++gg)
        {
            EdgeGlobal<U> & edge=graph.globals[gg];
            EdgeGlobalCopy<U> & edgecopy=this->globalcopys[gg];
            for(int a=0; a<tAsz-1; ++a)
                edgecopy.H1[a]=alloc_initial(edge.id, O.size(),0);
        }
        for(size_t w=0; w<=O.size(); ++w)
        {
            for(int ll=l; ll<l+graph.local_C[scc]; ++ll)
                CircleIterationEdge(w,ll,O,tAsz);
            for(int gg=g; gg<g+graph.global_C[scc]; ++gg)
                CircleIterationEdgeGlobal(w,gg,O,tAsz);
            for(int a=1; a<tAsz; ++a)
            {
                for(int ll=l; ll<l+graph.local_C[scc]; ++ll)
                {
                    EdgeLocal<U> & edge=graph.locals[ll];
                    EdgeLocalCopy<U> & edgecopy=this->localcopys[ll];
                    NodeCopy<U> & NTCP=this->nodecopys[edge.tail];
                    source_max(edgecopy.H1a[a-1][w],{&NTCP.tildeA[a-1][w]},{NTCP.tildeA[a-1][w].val+edge.T});
                    double gfwX=0;
                    if(edge.seq.size()>0)
                        gfwX=edge.uf[w]*(edge.seq.size()-1)+edge.vf[w];
                    source_max(edgecopy.H1b[a-1][w],{&edgecopy.H1a[a-1][w],&NTCP.tildeA[a-1][w]},{edgecopy.H1a[a-1][w].val+gfwX,NTCP.tildeA[a-1][w].val+edge.gfp[edge.seq.size()]+edge.T});
                    source_max(edgecopy.A[a][w],{&edgecopy.H1a[a-1][w],&edgecopy.H1b[a-1][w],&edgecopy.A[0][w]},{edgecopy.H1a[a-1][w].val+edge.gfm[0],edgecopy.H1b[a-1][w].val,edgecopy.A[0][w].val});
                }
                for(int gg=g; gg<g+graph.global_C[scc]; ++gg)
                {
                    EdgeGlobal<U> & edge=graph.globals[gg];
                    EdgeGlobalCopy<U> & edgecopy=this->globalcopys[gg];
                    NodeCopy<U> & NTCP=this->nodecopys[edge.tail];
                    source_max(edgecopy.H1[a-1][w],{&NTCP.tildeA[a-1][w]},{NTCP.tildeA[a-1][w].val+edge.T});
                    source_max(edgecopy.A[a][w],{&edgecopy.H1[a-1][w],&edgecopy.A[0][w]},{edgecopy.H1[a-1][w].val,edgecopy.A[0][w].val});
                }
                for(int i=0; i<graph.n_sz; ++i)
                    if(graph.nodes[i].scc==scc)
                    {
                        Dot<U> & tildeA=this->nodecopys[i].tildeA[a][w];
                        adrs.clear();
                        vals.clear();
                        if(this->nodecopys[i].barA)
                        {
                            adrs.push_back(&this->nodecopys[i].barA[w]);
                            vals.push_back(this->nodecopys[i].barA[w].val);
                        }
                        for(Dot<U>** ptr : this->nodecopys[i].ciptr)
                        {
                            adrs.push_back(&ptr[0][w]);
                            vals.push_back(ptr[0][w].val);
                        }
                        for(Dot<U>** ptr : this->nodecopys[i].siptr)
                        {
                            adrs.push_back(&ptr[a][w]);
                            vals.push_back(ptr[a][w].val);
                        }
                        source_max(tildeA,adrs.begin(),vals.begin(),vals.end());
                    }
            }
            if(w<O.size())
            {
                for(int ll=l; ll<l+graph.local_C[scc]; ++ll)
                {
                    EdgeLocal<U> & edge=graph.locals[ll];
                    EdgeLocalCopy<U> & edgecopy=this->localcopys[ll];
                    Node<U> & NH=graph.nodes[edge.head];
                    NodeCopy<U> & NHCP=this->nodecopys[edge.head];
                    source_max(edgecopy.B[w+1],{&edgecopy.B[w],&edgecopy.A[tAsz-1][w]},{edgecopy.B[w].val+NH.tue,edgecopy.A[tAsz-1][w].val+NH.tve});
                    if(edgecopy.B[w+1].val>NHCP.tildeA[0][w+1].val)
                    {
                        NHCP.tildeA[0][w+1].val=edgecopy.B[w+1].val;
                        NHCP.tildeA[0][w+1].s_sz=0;
                    }
                }
                for(int gg=g; gg<g+graph.global_C[scc]; ++gg)
                {
                    EdgeGlobal<U> & edge=graph.globals[gg];
                    EdgeGlobalCopy<U> & edgecopy=this->globalcopys[gg];
                    Node<U> & NH=graph.nodes[edge.head];
                    NodeCopy<U> & NHCP=this->nodecopys[edge.head];
                    source_max(edgecopy.B[w+1],{&edgecopy.B[w],&edgecopy.A[tAsz-1][w]},{edgecopy.B[w].val+NH.tue,edgecopy.A[tAsz-1][w].val+NH.tve});
                    if(edgecopy.B[w+1].val>NHCP.tildeA[0][w+1].val)
                    {
                        NHCP.tildeA[0][w+1].val=edgecopy.B[w+1].val;
                        NHCP.tildeA[0][w+1].s_sz=0;
                    }
                }
            }
        }
    }
    
    void CircleIterationEdge(int w, int l, std::string & O, int tAsz)
    {
        EdgeLocal<U> & edge=graph.locals[l];
        EdgeLocalCopy<U> & edgecopy=this->localcopys[l];
        NodeCopy<U> & NTCP=this->nodecopys[edge.tail];
        NodeCopy<U> & NHCP=this->nodecopys[edge.head];
        Dot<U> **E=edgecopy.E, **F=edgecopy.F, **G=edgecopy.G, *B=edgecopy.B, **C=edgecopy.C, **H2=edgecopy.H2, **A=edgecopy.A;
        double *ue=edge.ue, *ve=edge.ve, T=edge.T, *gfp=edge.gfp;
        std::vector<double> &uf=edge.uf, &vf=edge.vf;
        std::string & seq=edge.seq;
        if(w==0)
        {
            for(size_t s=0; s<=seq.size(); ++s)
            {
                E[s][w].val=-inf;
                F[s][w].val=-inf;
                G[s][w].val=-inf;
            }
            B[w].val=-inf;
            A[0][w].val=-inf;
        }
        else
        {
            C[0][w-1].val=-inf;
            source_max(H2[0][w-1],{&NTCP.tildeA[tAsz-1][w-1]},{NTCP.tildeA[tAsz-1][w-1].val+T});
            CircleInitial(w,E,F,G,H2[0][w-1],ve[0],ue[0]);
            for(size_t s=1; s<=seq.size(); ++s)
            {
                source_max(C[s][w-1],{&C[s-1][w-1],&H2[s-1][w-1]},{C[s-1][w-1].val+uf[w-1],H2[s-1][w-1].val+vf[w-1]});
                source_max(H2[s][w-1],{&C[s][w-1],&NTCP.tildeA[tAsz-1][w-1]},{C[s][w-1].val,NTCP.tildeA[tAsz-1][w-1].val+gfp[s]+T});
                CircleBody(seq,O,w,s,E,F,G,H2[s-1][w-1],H2[s][w-1],ve[s],ue[s],edge);
            }
            adrs.clear();
            vals.clear();
            adrs.push_back(&B[w]);
            vals.push_back(B[w].val);
            for(size_t s=0; s<=seq.size(); ++s)
            {
                adrs.push_back(&G[s][w]);
                vals.push_back(G[s][w].val);
            }
            source_max(A[0][w],adrs.begin(),vals.begin(),vals.end());
            if(A[0][w].val>=NHCP.tildeA[0][w].val)
            {
                adrs.clear();
                vals.clear();
                for(int j=0; j<NHCP.tildeA[0][w].s_sz; ++j)
                {
                    adrs.push_back(NHCP.tildeA[0][w].sources[j]);
                    vals.push_back(NHCP.tildeA[0][w].val);
                }
                adrs.push_back(&A[0][w]);
                vals.push_back(A[0][w].val);
                source_max(NHCP.tildeA[0][w],adrs.begin(),vals.begin(),vals.end());
            }
        }
    }
    
    void CircleInitial(int w, Dot<U>** E, Dot<U>** F, Dot<U>** G, Dot<U> & H2, double ve, double ue)
    {
        source_max(E[0][w],{&E[0][w-1],&H2},{E[0][w-1].val+ue,H2.val+ve});
        F[0][w].val=-inf;
        source_max(G[0][w],{&E[0][w]},{E[0][w].val});
    }
    
    template <typename Y>
    void CircleBody(std::string & seq, std::string & O, int w, int s, Dot<U>** E, Dot<U>** F, Dot<U>** G, Dot<U> & Hsm1, Dot<U> & Hs, double ve, double ue, Y & edge)
    {
        source_max(E[s][w],{&G[s][w-1],&E[s][w-1],&Hs},{G[s][w-1].val+ve,E[s][w-1].val+ue,Hs.val+ve});
        source_max(F[s][w],{&G[s-1][w],&F[s-1][w]},{G[s-1][w].val+edge.vf[w],F[s-1][w].val+edge.uf[w]});
        double gm=edge.gamma[int(seq[s-1])][int(O[w-1])];
        source_max(G[s][w],{&E[s][w],&F[s][w],&G[s-1][w-1],&Hsm1},{E[s][w].val,F[s][w].val,G[s-1][w-1].val+gm,Hsm1.val+gm});
    }
    
    void CircleIterationEdgeGlobal(int w, int g, std::string & O, int tAsz)
    {
        EdgeGlobal<U> & edge=graph.globals[g];
        EdgeGlobalCopy<U> & edgecopy=this->globalcopys[g];
        NodeCopy<U> & NTCP=this->nodecopys[edge.tail];
        Node<U> & NH=graph.nodes[edge.head];
        NodeCopy<U> & NHCP=this->nodecopys[edge.head];
        std::list<SNC<U>> &sncs=edgecopy.sncs;
        Dot<U> *B=edgecopy.B, **A=edgecopy.A;
        BroWheel<U> & Bro=edge.Bro;
        double T=edge.T, ue=edge.ue, ve=edge.ve, gamma_max=graph.gamma_max;
        std::vector<double> &uf=edge.uf, &vf=edge.vf;
        std::vector<std::vector<U>> &boxes=edgecopy.boxes;
        if(w==0)
        {
            sncs.emplace_back();
            sncs.back().E=sncs.back().F=sncs.back().G=-inf;
            sncs.back().sr1=0, sncs.back().sr2=Bro.size()-1;
            sncs.back().s=0;
            sncs.back().itd=sncs.back().itj=sncs.end();
            B[0].val=-inf;
            A[0][w].val=-inf;
        }
        else
        {
            double opt=B[w].val+(O.size()-w)*NH.tue;
            double H2=NTCP.tildeA[tAsz-1][w-1].val+T;
            auto it=sncs.begin();
            it->Gp=it->G;
            it->F=-inf;
            A[0][w].val=it->G=it->E=std::max(it->E+ue,H2+ve);
            boxes[w].clear();
            boxes[w].push_back(it->sr1);
            boxes[w].push_back(it->sr2);
            boxes[w].push_back(it->s);
            std::stack<typename std::list<SNC<U>>::iterator> to_add;
            InsertSNC(it,to_add,sncs,Bro);
            it=add_it(it,to_add,sncs.end());
            while(it!=sncs.end())
            {
                it->Gp=it->G;
                it->E=std::max(it->G+ve,it->E+ue);
                it->F=std::max(it->itp->G+vf[w],it->itp->F+uf[w]);
                double gm=edge.gamma[it->letter][int(O[w-1])];
                it->G=std::max(std::max(it->E,it->F),it->itp->Gp+gm);
                if(it->s==1)
                    it->G=std::max(it->G,H2+gm);
                if(it->G>=A[0][w].val)
                {
                    if(it->G>A[0][w].val)
                    {
                        A[0][w].val=it->G;
                        boxes[w].clear();
                        opt=std::max(opt,A[0][w].val+NH.ge[O.size()-w]);
                    }
                    boxes[w].push_back(it->sr1);
                    boxes[w].push_back(it->sr2);
                    boxes[w].push_back(it->s);
                }
                if(it->G<NTCP.tildeA[0][w].val+T)
                    it->G=it->F=-inf;
                if(it->E+ue<NTCP.tildeA[0][w].val+T+ve)
                    it->E=-inf;
                if(it->G+(O.size()-w)*gamma_max<opt)
                    it->G=it->F=it->E=-inf;
                InsertSNC(it,to_add,sncs,Bro);
                it=add_it(it,to_add,sncs.end());
            }
            if(size_t(w)==O.size())
                sncs.clear();
            else
                ShearTree(sncs.begin(),sncs.end());
            if(B[w].val>=A[0][w].val)
            {
                if(B[w].val>A[0][w].val)
                    boxes[w].clear();
                source_max(A[0][w],{&B[w]},{B[w].val});
            }
            else
                A[0][w].s_sz=0;
            if(A[0][w].val>=NHCP.tildeA[0][w].val)
            {
                adrs.clear();
                vals.clear();
                for(int j=0; j<NHCP.tildeA[0][w].s_sz; ++j)
                {
                    adrs.push_back(NHCP.tildeA[0][w].sources[j]);
                    vals.push_back(NHCP.tildeA[0][w].val);
                }
                adrs.push_back(&A[0][w]);
                vals.push_back(A[0][w].val);
                source_max(NHCP.tildeA[0][w],adrs.begin(),vals.begin(),vals.end());
            }
        }
    }
    
    void InsertSNC(typename std::list<SNC<U>>::iterator it, std::stack<typename std::list<SNC<U>>::iterator> & to_add, std::list<SNC<U>> & sncs, BroWheel<U> & Bro)
    {
        if(it->Gp!=-inf || it->G!=-inf)
        {
            if(it->itcs.size()==0)
            {
                auto itn=std::next(it);
                for(U letter=2; letter<7; ++letter)
                {
                    SNC<U> snc;
                    snc.sr1=it->sr1;
                    snc.sr2=it->sr2;
                    Bro.PreRange(snc.sr1, snc.sr2, letter);
                    if(snc.sr1<=snc.sr2)
                    {
                        snc.s=it->s+1;
                        snc.letter=letter;
                        snc.E=snc.F=snc.G=-inf;
                        snc.itp=it;
                        sncs.insert(itn,snc);
                        it->itcs.push_back(std::prev(itn));
                    }
                }
                if(it->itcs.size()==0)
                    it->itcs.push_back(sncs.end());
            }
            if(it->itcs[0]!=sncs.end())
                for(auto vit=it->itcs.rbegin(); vit!=it->itcs.rend(); ++vit)
                    to_add.push(*vit);
        }
    }
    
    typename std::list<SNC<U>>::iterator add_it(typename std::list<SNC<U>>::iterator it, std::stack<typename std::list<SNC<U>>::iterator> & to_add, typename std::list<SNC<U>>::iterator end_it)
    {
        if(!to_add.empty())
        {
            bool sg=(it->itd==end_it || to_add.top()->s>it->itd->s), se=(to_add.top()->s==it->itd->s), ll=(to_add.top()->letter<it->itd->letter), le=(to_add.top()->letter==it->itd->letter);
            if(sg || (se && (ll || le)))
            {
                if(sg || ll)
                {
                    to_add.top()->itd=it->itd;
                    it->itd=to_add.top();
                    to_add.top()->itj=it->itj;
                }
                it->itj=to_add.top();
                to_add.pop();
            }
        }
        return it->itj;
    }
    
    void ShearTree(typename std::list<SNC<U>>::iterator beg_it, typename std::list<SNC<U>>::iterator end_it)
    {
        auto it=beg_it;
        std::stack<typename std::list<SNC<U>>::iterator> down_stack;
        while(it->itd!=end_it)
        {
            down_stack.push(it->itd);
            auto itt=it->itj;
            while(itt!=end_it && itt->E==-inf && itt->G==-inf && itt->itp->G==-inf)
            {
                while(!down_stack.empty() && down_stack.top()->s>=itt->itd->s)
                    down_stack.pop();
                down_stack.push(itt->itd);
                itt->Gp=-inf;
                itt=itt->itj;
            }
            if(itt==end_it)
            {
                it->itd=end_it;
                it->itj=end_it;
                break;
            }
            else
            {
                auto ite=down_stack.top();
                do
                {
                    auto itop=down_stack.top();
                    down_stack.pop();
                    if(down_stack.empty())
                        it->itd=itop;
                    else
                        itop->itp->itd=itop;
                }while(!down_stack.empty());
                auto iti=it;
                while(iti!=ite)
                {
                    iti->itj=itt;
                    iti=iti->itd;
                }
                it=itt;
            }
        }
    }
    
    void update_cross(NodeCopy<U> & nodecopy, int W)
    {
        if(size_t(nodecopy.cid_now)==nodecopy.ciptr.size())
            for(int w=0; w<=W; ++w)
            {
                adrs.clear();
                vals.clear();
                if(nodecopy.barA)
                {
                    adrs.push_back(&nodecopy.barA[w]);
                    vals.push_back(nodecopy.barA[w].val);
                }
                for(Dot<U>** ptr : nodecopy.ciptr)
                {
                    adrs.push_back(&ptr[0][w]);
                    vals.push_back(ptr[0][w].val);
                }
                source_max(nodecopy.tildeA[0][w],adrs.begin(),vals.begin(),vals.end());
            }
    }
    
    void alloc_initial(Dot<U>** EFG, int eid, int S, int W, U* rob)
    {      
        for(int s=0; s<=S; ++s)
        {
            EFG[s]=heap_alloc<Dot<U>>(W+1);
            for(int w=0; w<=W; ++w)
            {
                EFG[s][w].eid=eid;
                EFG[s][w].s=s;
                EFG[s][w].w=w;
                EFG[s][w].rob=rob;
                EFG[s][w].visit=false;
            }
        }
    }
    
    Dot<U>* alloc_initial(int eid, int W)
    {
        Dot<U>* AB=heap_alloc<Dot<U>>(W+1);    
        for(int w=0; w<=W; ++w)
        {
            AB[w].eid=eid;
            AB[w].w=w;
            AB[w].rob=NULL;
            AB[w].visit=false;
        }
        return AB;
    }
    
    Dot<U>* alloc_initial(int eid, int W, int s)
    {
        Dot<U>* AB=alloc_initial(eid, W);    
        for(int w=0; w<=W; ++w)
            AB[w].s=s; // to distinguish AB from EFG, use s=-1 for A and s=-2 for B
        return AB;
    }
    
    void source_max(Dot<U> & dot, std::initializer_list<Dot<U>*> adrs, std::initializer_list<double> vals)
    {
        source_max(dot, adrs.begin(), vals.begin(), vals.end());
    }
    
    template <typename Ita, typename Itv>
    void source_max(Dot<U> & dot, Ita ita_b, Itv itv_b, Itv itv_e)
    {
        dot.val=-inf;
        dot.s_sz=0;
        for(Itv itv=itv_b; itv!=itv_e; ++itv)
        {
            if(*itv>dot.val)
                dot.val=*itv, dot.s_sz=1;
            else
                if(*itv==dot.val)
                    ++dot.s_sz;
        }
        dot.sources=heap_alloc<Dot<U>*>(dot.s_sz);
        int j=0;
        Ita ita=ita_b;
        for(Itv itv=itv_b; itv!=itv_e; ++itv,++ita)
            if(*itv==dot.val)
                dot.sources[j++]=*ita;
    }
};

#endif

