#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <cfloat>
#include <sstream>
#include "BroWheel.h"

constexpr static const double inf=std::numeric_limits<double>::infinity();

struct Edge
{
    int id;
    int tail;
    int head;
    bool global;
    std::list<std::string> files;
    bool visit;
    bool in_stack;
    int disc;
    int low;
    
    int id_tmp;
    int tail_tmp;
    int head_tmp;
    bool global_tmp;
    std::list<std::string> files_tmp;
    
    void temp(Edge & edge)
    {
        edge.id_tmp=id;
        edge.tail_tmp=tail;
        edge.head_tmp=head;
        edge.global_tmp=global;
        std::swap(edge.files_tmp,files);
    }
    
    void update()
    {
        id=id_tmp;
        tail=tail_tmp;
        head=head_tmp;
        global=global_tmp;
        std::swap(files,files_tmp);
    }
    
    void copy(Edge & edge)
    {
        id=edge.id;
        tail=edge.tail;
        head=edge.head;
        global=edge.global;
        std::swap(files,edge.files);
    }
    
    void swap(Edge & edge)
    {
        std::swap(id,edge.id);
        std::swap(tail,edge.tail);
        std::swap(head,edge.head);
        std::swap(global,edge.global);
        std::swap(files,edge.files);
    }
};

void getlineDoubles(std::vector<double> & container, std::ifstream & fin)
{
    std::string str_tmp;
    double tmp;
    std::getline(fin,str_tmp);
    std::istringstream iss(str_tmp);
    container.clear();
    iss >> tmp;
    while(!iss.fail())
    {
        container.push_back(tmp);
        iss >> tmp;
    }
}

template <typename T>
struct Dot // cannot use constructor because on Memory
{
    int eid; // eid=-1: Dot belongs node but not edge
    int s; // s=-1: Dot is A; s=-2: Dot is B.
    int w;
    double val;
    Dot<T>** sources;
    int s_sz;
    T* rob; // null if no r
    int id;
    bool visit; // visit must be false before backtrack
};

template <typename T>
struct SNC
{
    double E;
    double F;
    double G;
    double Gp;
    T letter;
    T sr1;
    T sr2;
    T s;
    typename std::list<SNC<T>>::iterator itp;
    std::vector<typename std::list<SNC<T>>::iterator> itcs;
    typename std::list<SNC<T>>::iterator itd;
};

template <typename T>
struct Node
{
    int tAsz=0;
    double tve;
    double tue;
    int scc;
    std::vector<double> ge;
};

template <typename U>
struct EdgeLocal
{
    int id;
    int tail;
    int head;
    std::string seq;
    double* ve=NULL;
    double* ue;
    double* gfp;
    double* gfm;
    double T;
    std::vector<double> vf;
    std::vector<double> uf;
    double gamma[7][7];
    
    void readin(std::string local_file, std::string score_file)
    {
        std::ifstream fin(local_file);
        fin >> seq;
        seq.shrink_to_fit();
        for(size_t i=0; i<seq.size(); ++i)
            seq[i]=(seq[i]/2+3) & 7;
        fin.close();
        
        fin.open(score_file);
        delete[] ve;
        ve=new double[4*seq.size()+4];
        ue=ve+seq.size()+1;
        gfp=ue+seq.size()+1;
        gfm=gfp+seq.size()+1;
        for(size_t i=0; i<=seq.size(); ++i)
            fin >> ve[i];
        for(size_t i=0; i<=seq.size(); ++i)
            fin >> ue[i];
        double tvfp, tufp, tvfm, tufm;
        fin >> tvfp >> tufp >> tvfm >> tufm >> T >> std::ws;
        gfp[0]=0, gfp[1]=tvfp, gfm[seq.size()]=0, gfm[seq.size()-1]=tvfm;
        for(size_t i=2,j=seq.size()-2; i<=seq.size(); ++i,--j)
            gfp[i]=gfp[i-1]+tufp, gfm[j]=gfm[j+1]+tufm;
        getlineDoubles(vf, fin);
        getlineDoubles(uf, fin);
        for(int i=2; i<7; ++i)
            for(int j=2; j<7; ++j)
                fin >> gamma[i][j];
        fin.close();
    }
    
    ~EdgeLocal()
    {
        delete[] ve;
    }
};

template <typename U>
struct EdgeGlobal
{
    int id;
    int tail;
    int head;
    BinWord<U> G;
    BroWheel<U> Bro;
    double ve;
    double ue;
    double T;
    std::vector<double> vf;
    std::vector<double> uf;
    double gamma[7][7];
    
    void readin_and_index(std::string index_file, std::string score_file)
    {
        std::ifstream fin(score_file);
        fin >> ve >> ue >> T >> std::ws;
        getlineDoubles(vf, fin);
        getlineDoubles(uf, fin);
        for(int i=2; i<7; ++i)
            for(int j=2; j<7; ++j)
                fin >> gamma[i][j];
        fin.close();
        Bro.reindex(G, index_file);
    }
};

template <typename T>
struct Graph
{
    int n_sz;
    Node<T>* nodes=NULL;
    int r_sz;
    int* roots=NULL;
    int t_sz;
    int* targets;
    int l_sz;
    EdgeLocal<T>* locals=NULL;
    int g_sz;
    EdgeGlobal<T>* globals=NULL;
    EdgeGlobal<T>** eid2globals=NULL;
    int scc_n;
    int* local_C=NULL;
    int* global_C;
    int* tAsz;
    bool* cross=NULL;
    
    Graph(std::string file)
    {
        std::ifstream fin(file);
        fin >> n_sz >> r_sz >> t_sz >> l_sz >> g_sz;
        delete[] nodes;
        nodes=new Node<T>[n_sz];
        for(int i=0; i<n_sz; ++i)
            fin >> nodes[i].tve >> nodes[i].tue;
        delete[] roots;
        roots=new int[r_sz+t_sz];
        targets=roots+r_sz;
        for(int i=0; i<r_sz; ++i)
            fin >> roots[i];
        for(int i=0; i<t_sz; ++i)
            fin >> targets[i];
        int e_sz=l_sz+g_sz;
        delete[] eid2globals;
        eid2globals=new EdgeGlobal<T>*[e_sz];
        delete[] cross;
        cross=new bool[e_sz];
        Edge* edges=new Edge[e_sz];
        for(int i=0; i<e_sz; ++i)
        {
            eid2globals[i]=NULL;
            edges[i].id=i;
            fin >> edges[i].tail >> edges[i].head;
            edges[i].global=(i>=l_sz);
            edges[i].files.resize(2);
            fin >> edges[i].files.front() >> edges[i].files.back();
        }
        ShearSort(edges);
        delete[] locals;
        locals=new EdgeLocal<T>[l_sz];
        delete[] globals;
        globals=new EdgeGlobal<T>[g_sz];
        for(int i=0, e=0, l=0, g=0; i<scc_n; ++i)
        {
            for(int j=0; j<local_C[i]; ++j,++e,++l)
            {
                locals[l].id=edges[e].id;
                locals[l].tail=edges[e].tail;
                locals[l].head=edges[e].head;
                locals[l].readin(edges[e].files.front(),edges[e].files.back());
                
                if(nodes[locals[l].tail].tAsz<tAsz[i])
                    nodes[locals[l].tail].tAsz=tAsz[i];
                if(local_C[i]==1 && global_C[i]==0 && locals[l].tail!=locals[l].head)
                {
                    cross[locals[l].id]=true;
                }
                else
                {
                    cross[locals[l].id]=false;
                    nodes[locals[l].head].scc=i;
                }
            }
            for(int j=0; j<global_C[i]; ++j,++e,++g)
            {
                globals[g].id=edges[e].id;
                globals[g].tail=edges[e].tail;
                globals[g].head=edges[e].head;
                globals[g].readin_and_index(edges[e].files.front(),edges[e].files.back());
                
                if(nodes[globals[g].tail].tAsz<tAsz[i])
                    nodes[globals[g].tail].tAsz=tAsz[i];
                if(local_C[i]==0 && global_C[i]==1 && globals[g].tail!=globals[g].head)
                {
                    cross[globals[g].id]=true;
                }
                else
                {
                    cross[globals[g].id]=false;
                    nodes[globals[g].head].scc=i;
                }
                eid2globals[globals[g].id]=&globals[g];
            }
        }
        delete[] edges;
        for(int i=0; i<t_sz; ++i)
            if(nodes[targets[i]].tAsz<1)
                nodes[targets[i]].tAsz=1;
        for(int i=0; i<n_sz; ++i)
            if(nodes[i].tAsz>0)
            {
                nodes[i].ge.push_back(0);
                nodes[i].ge.push_back(nodes[i].tve);
            }
    }
    
    ~Graph()
    {
        delete[] nodes;
        delete[] roots;
        delete[] locals;
        delete[] globals;
        delete[] eid2globals;
        delete[] local_C;
        delete[] cross;
    }
    
    void ShearSort(Edge* edges)
    {
        RemoveEdge(edges, roots, r_sz);
        int e_sz=l_sz+g_sz;
        for(int i=0; i<e_sz; ++i)
            std::swap(edges[i].tail,edges[i].head);
        RemoveEdge(edges, targets, t_sz);
        e_sz=l_sz+g_sz;
        std::unordered_set<int> nzout;
        for(int i=0; i<e_sz; ++i)
        {
            std::swap(edges[i].tail,edges[i].head);
            nzout.insert(edges[i].tail);
        }
        int new_rt_sz=0;
        for(int i=0; i<r_sz; ++i)
            if(nzout.find(roots[i])!=nzout.end())
            {
                if(i!=new_rt_sz)
                    roots[new_rt_sz]=roots[i];
                ++new_rt_sz;
            }
        r_sz=new_rt_sz;
        
        std::list<int> C;
        SCCTS(edges, e_sz, C);
        scc_n=C.size();
        std::unordered_map<int,int> indegree;
        for(int i=0; i<e_sz; ++i)
            AccumInsert(indegree,edges[i].head);
        delete[] local_C;
        local_C=new int[3*scc_n];
        global_C=local_C+scc_n;
        tAsz=global_C+scc_n;
        int it=0;
        for(int e_szg : C)
        {
            int* heads=new int[e_szg];
            for(int i=0; i<e_szg; ++i)
                heads[i]=edges[i].head;
            std::sort(heads,heads+e_szg);
            tAsz[it]=1;
            for(int i=1; i<e_szg; ++i)
                if(heads[i]!=heads[i-1])
                    ++tAsz[it];
            delete[] heads;
            int s=0, e=e_szg-1;
            while(true)
            {
                while(s<=e && !edges[s].global)
                    ++s;
                while(s<=e && edges[e].global)
                    --e;
                if(s<e)
                    edges[s].swap(edges[e]);
                else
                    break;
            }
            local_C[it]=s;
            global_C[it]=e_szg-s;
            edges+=s;
            e_szg-=s;
            GeneralizedTS(edges, e_szg, indegree);
            edges+=e_szg;
            ++it;
        }
    }

    void GeneralizedTS(Edge* edges, int e_sz, std::unordered_map<int,int> & indegree)
    {
        if(e_sz<2)
            return;
        std::list<int> C;
        SCCTS(edges, e_sz, C);
        for(int e_szg : C)
        {
            if(e_szg==1)
            {
                ++edges;
                continue;
            }
            std::unordered_map<int,int> indegreeSCC,outdegreeSCC;
            for(int i=0; i<e_szg; ++i)
            {
                AccumInsert(indegreeSCC,edges[i].head);
                AccumInsert(outdegreeSCC,edges[i].tail);
            }
            double MaxRatio=-std::numeric_limits<double>::infinity();
            int e;
            for(int i=0; i<e_szg; ++i)
            {
                int tmp=indegreeSCC[edges[i].tail];
                double NowRatio=double(outdegreeSCC[edges[i].head])/double(tmp);
                if(indegree[edges[i].tail]>tmp && NowRatio>MaxRatio)
                {
                    MaxRatio=NowRatio;
                    e=i;
                }
            }
            edges[0].swap(edges[e]);
            GeneralizedTS(++edges, --e_szg, indegree);
            edges+=e_szg;
        }
    }
    
    void AccumInsert(std::unordered_map<int,int> & degree, int key)
    {
        auto iter_bool=degree.insert(std::make_pair(key,1));
        if(!(iter_bool.second))
            ++(iter_bool.first->second);
    }

    void SCCTS(Edge* edges, int e_sz, std::list<int> & C)
    {
        for(int i=0; i<e_sz; ++i)
            edges[i].visit=false;
        std::unordered_map<int,std::list<int>> lists;
        for(int i=0; i<e_sz; ++i)
            lists[edges[i].tail].push_back(i);
        std::stack<int> S;
        for(int i=0; i<e_sz; ++i)
            edges[i].in_stack=false;
        Edge* edge_E=edges+e_sz;
        int id=0;
        for(int i=0; i<e_sz; ++i)
            if(!edges[i].visit)
                DeepFirst(i, edges, e_sz, lists, S, edge_E, C, id);
        for(int i=0; i<e_sz; ++i)
            edges[i].update();
    }
    
    void RemoveEdge(Edge* edges, int* rt, int & rt_sz)
    {
        int e_sz=l_sz+g_sz;
        for(int i=0; i<e_sz; ++i)
            edges[i].visit=false;
        std::unordered_map<int,std::list<int>> lists;
        for(int i=0; i<e_sz; ++i)
            lists[edges[i].tail].push_back(i);
        int new_rt_sz=0;
        for(int i=0; i<rt_sz; ++i)
            if(lists.find(rt[i])!=lists.end())
            {
                if(i!=new_rt_sz)
                    rt[new_rt_sz]=rt[i];
                ++new_rt_sz;
            }
        rt_sz=new_rt_sz;
        for(int i=0; i<rt_sz; ++i)
            for(int j : lists[rt[i]])
                if(!edges[j].visit)
                    DeepFirst(j, edges, e_sz, lists);
        
        int new_e_sz=0;
        for(int i=0; i<e_sz; ++i)
        {
            if(!edges[i].visit)
            {
                if(edges[i].global)
                    --g_sz;
                else
                    --l_sz;
            }
            else
            {
                if(i!=new_e_sz)
                    edges[new_e_sz].copy(edges[i]);
                ++new_e_sz;
            }
        }
        e_sz=new_e_sz;
    }

    void DeepFirst(int now, Edge* edges, int e_sz, std::unordered_map<int,std::list<int>> & lists, std::stack<int> & S, Edge* & edge_E, std::list<int> & C, int & id)
    {
        edges[now].visit=true;
        edges[now].disc=edges[now].low=id++;
        S.push(now);
        edges[now].in_stack=true;
        for(int i : lists[edges[now].head])
        {
            if(!edges[i].visit)
            {
                DeepFirst(i, edges, e_sz, lists, S, edge_E, C, id);
                edges[now].low=std::min(edges[now].low,edges[i].low);
            }
            else
            {
                if(edges[i].in_stack)
                    edges[now].low=std::min(edges[now].low,edges[i].disc);
            }
        }
        if(edges[now].disc==edges[now].low)
        {
            int tmp;
            C.push_front(0);
            do
            {
                tmp=S.top();
                S.pop();
                edges[tmp].in_stack=false;
                edges[tmp].temp(*(--edge_E));
                ++C.front();
            }while(tmp!=now);
        }
    }
    
    void DeepFirst(int now, Edge* edges, int e_sz, std::unordered_map<int,std::list<int>> & lists)
    {
        edges[now].visit=true;
        for(int i : lists[edges[now].head])
            if(!edges[i].visit)
                DeepFirst(i, edges, e_sz, lists);
    }
};

template <typename T>
struct NodeCopy
{
    Dot<T>** tildeA=NULL;
    std::vector<Dot<T>**> ciptr;
    int cid_now=0;
    std::vector<Dot<T>**> siptr;
    Dot<T>* barA=NULL;
    
    ~NodeCopy()
    {
        delete[] tildeA;
    }
};

template <typename U>
struct EdgeLocalCopy
{
    Dot<U>** A=NULL;
    Dot<U>* B;
    Dot<U>** E;
    Dot<U>** F;
    Dot<U>** G;
    Dot<U>** H1a=NULL;
    Dot<U>** H1b;
    Dot<U>** H2;
    Dot<U>** C;
    
    ~EdgeLocalCopy()
    {
        delete[] A;
        delete[] H1a;
    }
};

template <typename U>
struct EdgeGlobalCopy
{
    std::vector<std::vector<U>> boxes;
    Dot<U>** A=NULL;
    Dot<U>* B;
    Dot<U>** H1=NULL;
    std::list<SNC<U>> sncs;
    
    ~EdgeGlobalCopy()
    {
        delete[] A;
        delete[] H1;
    }
};

template <typename T>
struct GraphCopy
{
    NodeCopy<T>* nodecopys=NULL;
    EdgeLocalCopy<T>* localcopys=NULL;
    EdgeGlobalCopy<T>* globalcopys=NULL;
    
    ~GraphCopy()
    {
        delete[] nodecopys;
        delete[] localcopys;
        delete[] globalcopys;
    }
    
    void GetCopy(Graph<T> & graph)
    {
        delete[] nodecopys;
        nodecopys=new NodeCopy<T>[graph.n_sz];
        delete[] localcopys;
        localcopys=new EdgeLocalCopy<T>[graph.l_sz];
        delete[] globalcopys;
        globalcopys=new EdgeGlobalCopy<T>[graph.g_sz];
        for(int i=0, e=0, l=0, g=0; i<graph.scc_n; ++i)
        {
            for(int j=0; j<graph.local_C[i]; ++j,++e,++l)
            {
                int EFGsz=graph.locals[l].seq.size()+1;
                localcopys[l].A=new Dot<T>*[graph.tAsz[i]+3*EFGsz];
                localcopys[l].E=localcopys[l].A+graph.tAsz[i];
                localcopys[l].F=localcopys[l].E+EFGsz;
                localcopys[l].G=localcopys[l].F+EFGsz;
                if(graph.local_C[i]==1 && graph.global_C[i]==0 && graph.locals[l].tail!=graph.locals[l].head)
                {
                    nodecopys[graph.locals[l].head].ciptr.push_back(localcopys[l].A);
                }
                else
                {
                    nodecopys[graph.locals[l].head].siptr.push_back(localcopys[l].A);
                    localcopys[l].H1a=new Dot<T>*[2*EFGsz+2*(graph.tAsz[i]-1)];
                    localcopys[l].H1b=localcopys[l].H1a+graph.tAsz[i]-1;
                    localcopys[l].H2=localcopys[l].H1b+graph.tAsz[i]-1;
                    localcopys[l].C=localcopys[l].H2+EFGsz;
                }
            }
            for(int j=0; j<graph.global_C[i]; ++j,++e,++g)
            {
                globalcopys[g].A=new Dot<T>*[graph.tAsz[i]];
                if(graph.local_C[i]==0 && graph.global_C[i]==1 && graph.globals[g].tail!=graph.globals[g].head)
                {
                    nodecopys[graph.globals[g].head].ciptr.push_back(globalcopys[g].A);
                }
                else
                {
                    nodecopys[graph.globals[g].head].siptr.push_back(globalcopys[g].A);
                    globalcopys[g].H1=new Dot<T>*[graph.tAsz[i]-1];
                }
            }
        }
        for(int i=0; i<graph.n_sz; ++i)
            if(graph.nodes[i].tAsz>0)
                nodecopys[i].tildeA=new Dot<T>*[graph.nodes[i].tAsz];
    }
};

#endif
