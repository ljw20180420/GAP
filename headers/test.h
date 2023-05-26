#ifndef TEST_H
#define TEST_H
#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "Align.h"

bool connect(int n_sz, std::set<int> & rt, std::multiset<std::pair<int,int>> & locals, std::multiset<std::pair<int,int>> & globals, bool root)
{
    std::vector<bool> visit(n_sz,false);
    for(int id : rt)
        visit[id]=true;
    int now_visit=rt.size(), last_visit;
    do
    {
        last_visit=now_visit;
        for(auto & edges : {locals,globals})
            for(auto & pair : edges)
            {
                int t=pair.first, h=pair.second;
                if(root)
                    visit[h]=(visit[h] || visit[t]);
                else
                    visit[t]=(visit[t] || visit[h]);
            }
        now_visit=0;
        for(int i=0; i<n_sz; ++i)
            if(visit[i])
                ++now_visit;
    }while(now_visit!=last_visit);
    return (now_visit==n_sz);
}

void plot_gamma(double mat, double mis, std::ofstream & fout)
{
    for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
        {
            char dil;
            double val;
            if(i>0 && i==j)
                val=mat;
            else
                val=mis;
            if(j<4)
                dil='\t';
            else
                dil='\n';
            fout << val << dil;
        }
    fout.close();
}

std::string random_mut(std::string & sstr, double indel_rate, double mut_rate, int head_in, int tail_in)
{
    std::string int2str("ATCG");
    std::string mut;
    size_t row=0;
    while(row<sstr.size())
    {
        double rv=double(rand())/RAND_MAX;
        if(rv<indel_rate)
            mut.push_back(int2str[rand()%4]);
        else if(rv<2*indel_rate)
            ++row;
        else
        {
            rv=double(rand())/RAND_MAX;
            if(rv<mut_rate)
                mut.push_back(int2str[rand()%4]);
            else
                mut.push_back(sstr[row++]);
        }
    }
    while(row<=sstr.size())
    {
        double rv=double(rand())/RAND_MAX;
        if(rv<indel_rate)
            mut.push_back(int2str[rand()%4]);
        else
            ++row;
    }
    std::string head_str, tail_str;
    head_in=rand()%head_in;
    tail_in=rand()%tail_in;
    for(int i=0; i<head_in; ++i)
        head_str.push_back(int2str[rand()%4]);
    for(int i=0; i<tail_in; ++i)
        tail_str.push_back(int2str[rand()%4]);
    return head_str+mut+tail_str;
}

void random_seq(int root, std::set<int> & targets, std::vector<std::pair<int,int>> & locals, std::vector<std::pair<int,int>> & globals, std::vector<std::string> & locals_seq, std::vector<std::string> & globals_seq, std::string & seq, std::ofstream & truth, int aseqlb, int asequb, double indel_rate, double mut_rate, int head_in, int tail_in, double apro)
{
    if(double(rand())/RAND_MAX<apro && targets.find(root)!=targets.end())
        return;
    else
    {
        std::vector<std::pair<int,int>> id_heads;
        size_t id=0;
        for(auto & edges : {locals,globals})
            for(auto & pair : edges)
            {
                if(pair.first==root)
                    id_heads.emplace_back(id,pair.second);
                ++id;
            }
        if(id_heads.empty())
            return;
        int rv=rand()%id_heads.size();
        id=id_heads[rv].first;
        int head=id_heads[rv].second;
        std::string* seq_ptr;
        if(id<locals_seq.size())
            seq_ptr=&locals_seq[id];
        else
            seq_ptr=&globals_seq[id-locals_seq.size()];
        size_t s,e;
        do
        {
            s=rand()%(seq_ptr->size()-aseqlb+1);
            e=s+aseqlb+rand()%(asequb-aseqlb);
        }while(e>seq_ptr->size());
        truth << id << '\t' << s << '\t' << e << '\n';
        std::string sstr=seq_ptr->substr(s,e-s);
        seq+=random_mut(sstr,indel_rate,mut_rate,head_in,tail_in);
        random_seq(head,targets,locals,globals,locals_seq,globals_seq,seq,truth,aseqlb,asequb,indel_rate,mut_rate,head_in,tail_in,apro);
    }
}

void random_DG(int n_sz, int r_sz, double gpro, double rpro, double tve, double tue, std::string graph_file, std::string local_file, std::string local_file_score, std::string index_file, std::string global_file, std::string global_file_score, int lseqlb, int lsequb, int gseqlb, int gsequb, std::string index_information, double ve, double ue, double vf, double uf, double T, double tvf, double tuf, double mat, double mis, int aseqlb, int asequb, int seq_num, std::string read_file, std::string truth_file, double indel_rate, double mut_rate, int head_in, int tail_in, bool acyclic, double apro)
{
    std::string int2str("ATCG");
    srand(1);
    std::multiset<std::pair<int,int>> locals;
    std::multiset<std::pair<int,int>> globals;
    std::set<int> roots;
    roots.insert(0);
    while(int(roots.size())<r_sz)
        roots.insert(rand()%(n_sz-1));
    std::set<int> targets;
    targets.insert(n_sz-1);
    bool flag;
    do
    {
        int head, tail;
        do
        {
            head=rand()%n_sz;
            tail=rand()%n_sz;
        }while(acyclic && head<=tail);
        std::pair<int,int> th=std::make_pair(tail,head);
        if((locals.count(th)==0 && globals.count(th)==0) || double(rand())/RAND_MAX<rpro)
        {
            if(double(rand())/RAND_MAX<gpro)
                globals.insert(th);
            else
                locals.insert(th);
        }
        else
            continue;
        flag=true;
        for(int i=0; i<n_sz; ++i)
            flag=flag && connect(n_sz, roots, locals, globals, true) && connect(n_sz, targets, locals, globals, false); 
    }while(!flag);
    
    std::ofstream fout(graph_file);
    fout << n_sz << '\n' << r_sz << '\n' << 1 << '\n' << locals.size() << '\n' << globals.size() << '\n';
    for(int i=0; i<n_sz; ++i)
        fout << tve << '\t' << tue << '\n';
    for(int i : roots)
        fout << i << '\n';
    fout << n_sz-1 << '\n';
    std::vector<std::string> locals_seq;
    int fn=0;
    for(auto & pair : locals)
    {
        fout << pair.first << '\t' << pair.second << '\t' << local_file+std::to_string(fn) << '\t' << local_file_score+std::to_string(fn) << '\n';
        std::ofstream lout(local_file+std::to_string(fn));
        int seqlen=lseqlb+rand()%(lsequb-lseqlb);
        locals_seq.emplace_back();
        for(int j=0; j<seqlen; ++j)
            locals_seq.back()+=int2str[rand()%4];
        lout << locals_seq.back() << '\n';
        lout.close();
        lout.open(local_file_score+std::to_string(fn));
        for(double xe : {ve,ue})
            for(int i=0; i<=seqlen; ++i)
            {
                lout << xe;
                if(i<seqlen)
                    lout << '\t';
                else
                    lout << '\n';
            }
        lout << tvf << '\t' << tuf << '\t' << tvf << '\t' << tuf << '\t' << T << '\n' << vf << '\n' << uf << '\n';
        plot_gamma(mat, mis, lout);
        lout.close();
        ++fn;
    }
    std::ofstream fout2(index_information);
    fout2 << globals.size() << '\n';
    std::vector<std::string> globals_seq;
    fn=0;
    for(auto & pair : globals)
    {
        fout << pair.first << '\t' << pair.second << '\t' << index_file+std::to_string(fn) << '\t' << global_file_score << '\n';
        fout2 << "1\n" << global_file+std::to_string(fn) << '\n' << index_file+std::to_string(fn) << '\n';
        std::ofstream gout(global_file+std::to_string(fn));
        globals_seq.emplace_back();
        int seqlen=gseqlb+rand()%(gsequb-gseqlb);
        for(int j=0; j<seqlen; ++j)
            globals_seq.back()+=int2str[rand()%4];
        gout << globals_seq.back() << '\n';
        gout.close();
        ++fn;
    }
    fout.close();
    fout2.close();
    fout.open(global_file_score);
    fout << ve << '\t' << ue << '\t' << T << '\n' << vf << '\n' << uf << '\n';
    plot_gamma(mat, mis, fout);
    fout.close();
    
    std::vector<std::pair<int,int>> locals_copy(locals.begin(),locals.end());
    std::vector<std::pair<int,int>> globals_copy(globals.begin(),globals.end());
    fout.open(read_file);
    fout2.open(truth_file);
    std::string seq;
    for(int sn=0; sn<seq_num; ++sn)
    {
        fout2 << sn << '\n';
        seq.clear();
        auto it=roots.begin();
        std::advance(it,rand()%roots.size());
        random_seq(*it,targets,locals_copy,globals_copy,locals_seq,globals_seq,seq,fout2,aseqlb,asequb,indel_rate,mut_rate,head_in,tail_in,apro);
        fout << seq << '\n';
    }
    fout.close();
    fout2.close();
}

void to_gv(std::string graph_file)
{
    std::ifstream fin(graph_file);
    std::ofstream fout(graph_file+".gv");
    fout << "digraph align_graph\n{\n";
    fout << "\tnode[shape=circle]\n";
    int n_sz, r_sz, t_sz, l_sz, g_sz;
    fin >> n_sz >> r_sz >> t_sz >> l_sz >> g_sz;
    std::string st;
    for(int i=0; i<n_sz*2; ++i)
        fin >> st;
    int node;
    for(int i=0; i<r_sz; ++i)
    {
        fin >> node;
        fout << '\t' << node << "[shape=square]\n";
    }
    for(int i=0; i<t_sz; ++i)
    {
        fin >> node;
        fout << '\t' << node << "[shape=diamond]\n";
    }
    for(int i=0,tail,head; i<l_sz+g_sz; ++i)
    {
        fin >> tail >> head >> st >> st;
        fout << '\t' << tail << "->" << head;
        if(i<l_sz)
            fout << "[color=black,label=\""+std::to_string(i)+"\"]\n";
        else
            fout << "[color=red,label=\""+std::to_string(i)+"\"]\n";
    }
    fout << "}\n";
    fin.close();
    fout.close();
    system(("dot -Tpdf "+graph_file+".gv > "+graph_file+".pdf").c_str());
}

#endif
