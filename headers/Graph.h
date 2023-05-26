#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <cfloat>
#include <list>
#include "BroWheel.h"
#include "Memory.h"

constexpr static const double inf = std::numeric_limits<double>::infinity();

struct Dot // cannot use constructor because on Memory
{
    int n; // n determines the type of Dot
    size_t s;
    int w;
    double val;
    Dot **sources; // apply on memory
    int s_sz;      // s_sz<0 means visited
    int lambda;
    size_t id;
};

struct Node;

struct Edge
{
    const static int sigma = 6;

    std::string name;
    double gamma[sigma + 1][sigma + 1];
    double ve;
    double ue;
    double vf;
    double uf;
    double T;
    int n;
    Node *tail;
    Node *head;

    Edge(std::string name_, double gamma_[7][7], double ve_, double ue_, double vf_, double uf_, double T_, int n_)
    {
        name = name_;
        for (int i = 0; i <= sigma; ++i)
            for (int j = 0; j <= sigma; ++j)
                gamma[i][j] = gamma_[i][j];
        ve = ve_;
        ue = ue_;
        vf = vf_;
        uf = uf_;
        T = T_;
        n = n_;
    }
};

struct NameSeq
{
    std::string name;
    std::string seq;

    void readin(std::string file)
    {
        std::ifstream fin(file);
        std::getline(std::getline(fin, name), seq);
    }
};

struct EdgeLocal : Edge
{
    NameSeq &nameseq;
    double vfp;
    double ufp;
    double vfm;
    double ufm;

    EdgeLocal(std::string name_, double gamma_[7][7], double ve_, double ue_, double vf_, double uf_, double T_, int n_, NameSeq &nameseq_, double vfp_, double ufp_, double vfm_, double ufm_)
        : Edge(name_, gamma_, ve_, ue_, vf_, uf_, T_, n_), nameseq(nameseq_)
    {
        vfp = vfp_;
        ufp = ufp_;
        vfm = vfm_;
        ufm = ufm_;
    }
};

struct EdgeLocalCross : EdgeLocal
{
    Dot **E;
    Dot **F;
    Dot **G;

    EdgeLocalCross(EdgeLocal &edge)
    : EdgeLocal(edge.name, edge.gamma, edge.ve, edge.ue, edge.vf, edge.uf, edge.T, edge.n, edge.nameseq, edge.vfp, edge.ufp, edge.vfm, edge.ufm)
    {
        tail=edge.tail;
        head=edge.head;
    }
};

struct EdgeLocalCircuit : EdgeLocal
{
    Dot *D;
    Dot **E;
    Dot **F0;
    Dot **G0;
    Dot **G;
    Dot **D0;
    Dot **DX;

    EdgeLocalCircuit(EdgeLocal &edge)
    : EdgeLocal(edge.name, edge.gamma, edge.ve, edge.ue, edge.vf, edge.uf, edge.T, edge.n, edge.nameseq, edge.vfp, edge.ufp, edge.vfm, edge.ufm)
    {
        tail=edge.tail;
        head=edge.head;
    }
};

struct TrackNode
{
    TrackNode* itp;
    TrackNode* itcs[5];
    Dot *E;
    Dot *F;
    Dot *G;
    int tau = 0;
    int lambda;

    TrackNode(Memory *memory_, TrackNode* itp_, TrackNode* itc_, int n, size_t s, int W, int lambda_)
    {
        itp = itp_;
        for (int i = 0; i < 5; ++i)
            itcs[i] = itc_;
        alloc_initial(memory_, E, n, s, W);
        alloc_initial(memory_, F, n, s, W);
        alloc_initial(memory_, G, n, s, W);
        lambda=lambda_;
        for (int w=0; w<=W; ++w)
        {
            E[w].lambda=lambda;
            F[w].lambda=lambda;
            G[w].lambda=lambda;
        }
    }

    void alloc_initial(Memory *memory_, Dot *&ptr, int n, size_t s, int W)
    {
        ptr = memory_->heap_alloc<Dot>(W+1);
        for (int w = 0; w <= W; ++w)
        {
            ptr[w].n = n;
            ptr[w].s = s;
            ptr[w].w = w;
        }
    }
};

struct EdgeGlobal : Edge
{
    BroWheel &browheel;
    std::deque<std::deque<std::pair<int64_t, int>>> Cdelta;
    std::deque<TrackNode> tracknodes;

    EdgeGlobal(std::string name_, double gamma_[7][7], double ve_, double ue_, double vf_, double uf_, double T_, int n_, BroWheel &browheel_)
        : Edge(name_, gamma_, ve_, ue_, vf_, uf_, T_, n_), browheel(browheel_)
    {
    }
};

struct EdgeGlobalCross : EdgeGlobal
{
    EdgeGlobalCross(EdgeGlobal &edge)
    : EdgeGlobal(edge.name, edge.gamma, edge.ve, edge.ue, edge.vf, edge.uf, edge.T, edge.n, edge.browheel)
    {
        tail=edge.tail;
        head=edge.head;
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
        : E(E_), F0(F0_), G0(G0_), hatG(hatG_), c(c_), idl(idl_), lambda(lambda_), sr1(sr1_), sr2(sr2_), cid(cid_), itp(itp_), jump(jump_)
    {
    }
};

struct EdgeGlobalCircuit : EdgeGlobal
{
    std::deque<SNC> sncs;

    Dot *G;
    Dot **D0;

    EdgeGlobalCircuit(EdgeGlobal &edge)
    : EdgeGlobal(edge.name, edge.gamma, edge.ve, edge.ue, edge.vf, edge.uf, edge.T, edge.n, edge.browheel)
    {
        tail=edge.tail;
        head=edge.head;
    }
};

struct Node
{
    std::string name;
    double ve;
    double ue;
    std::vector<EdgeLocalCross *> out_local_crosses;
    std::vector<EdgeLocalCircuit *> out_local_circuits;
    std::vector<EdgeGlobalCross *> out_global_crosses;
    std::vector<EdgeGlobalCircuit *> out_global_circuits;
    std::vector<EdgeLocalCross *> in_local_crosses;
    std::vector<EdgeLocalCircuit *> in_local_circuits;
    std::vector<EdgeGlobalCross *> in_global_crosses;
    std::vector<EdgeGlobalCircuit *> in_global_circuits;
    Dot Abar;
    std::vector<std::deque<Dot>> A;
    std::deque<Dot> B;
    std::deque<std::deque<Dot *>> AdeltaDot;
    std::deque<std::deque<EdgeGlobalCross *>> AdeltaCross;
    std::deque<std::deque<EdgeGlobalCircuit *>> AdeltaCircuit;

    Node(std::string name_, double ve_, double ue_)
    {
        name = name_;
        ve = ve_;
        ue = ue_;
    }
};

struct SCC
{
    std::vector<Node *> nodes;
    std::vector<EdgeLocalCross *> local_crosses;
    std::vector<EdgeLocalCircuit *> local_circuits;
    std::vector<EdgeGlobalCross *> global_crosses;
    std::vector<EdgeGlobalCircuit *> global_circuits;
};

struct Graph
{
    const static int sigma = 6;

    std::vector<Node> nodes;
    std::vector<Node *> roots;
    std::vector<Node *> targets;
    std::deque<EdgeLocalCross> local_crosses;
    std::deque<EdgeLocalCircuit> local_circuits;
    std::deque<EdgeGlobalCross> global_crosses;
    std::deque<EdgeGlobalCircuit> global_circuits;

    std::deque<SCC> sccs;

    Graph(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel)
    {
        for (int i = 1; i < argc; ++i)
            if (!strcmp(argv[i], "--nodes"))
            {
                while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                {
                    nodes.emplace_back(argv[i], str2double(argv[i + 1]), str2double(argv[i + 2]));
                    i += 2;
                }
                --i;
            }
        std::list<EdgeLocal> locals;
        std::list<EdgeGlobal> globals;
        for (int i = 1, n = 0; i < argc; ++i)
            if (!strcmp(argv[i], "--roots"))
            {
                while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                    for (auto &node : nodes)
                        if (node.name == argv[i])
                        {
                            roots.push_back(&node);
                            break;
                        }
                --i;
            }
            else if (!strcmp(argv[i], "--targets"))
            {
                while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                    for (auto &node : nodes)
                        if (node.name == argv[i])
                        {
                            targets.push_back(&node);
                            break;
                        }
                --i;
            }
            else if (!strcmp(argv[i], "--locals"))
            {
                while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                {
                    double gamma_[7][7]={
                        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
                        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
                        {-inf, -inf, -3, -3, -3, -3, -3},
                        {-inf, -inf, -3, 1, -3, -3, -3},
                        {-inf, -inf, -3, -3, 1, -3, -3},
                        {-inf, -inf, -3, -3, -3, 1, -3},
                        {-inf, -inf, -3, -3, -3, -3, 1}
                    };
                    if (!strcmp(argv[i],"default_gamma"))
                        ++i;
                    else
                        for (int a=0; a<7; ++a)
                            for (int b=0; b<7; ++b)
                                gamma_[a][b]=str2double(argv[i++]);
                    locals.emplace_back(argv[i], gamma_, str2double(argv[i + 1]), str2double(argv[i + 2]), str2double(argv[i + 3]), str2double(argv[i + 4]), str2double(argv[i + 5]), n++, file2seq[argv[i]], str2double(argv[i + 6]), str2double(argv[i + 7]), str2double(argv[i + 8]), str2double(argv[i + 9]));
                    i += 9;
                    for (auto &node : nodes)
                    {
                        if (node.name == argv[i + 1])
                            locals.back().tail = &node;
                        if (node.name == argv[i + 2])
                            locals.back().head = &node;
                    }
                    i += 2;
                }
                --i;
            }
            else if (!strcmp(argv[i], "--globals"))
            {
                while (++i < argc && (strlen(argv[i]) < 2 || argv[i][0] != '-' || argv[i][1] != '-'))
                {
                    double gamma_[7][7]={
                        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
                        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
                        {-inf, -inf, -3, -3, -3, -3, -3},
                        {-inf, -inf, -3, 1, -3, -3, -3},
                        {-inf, -inf, -3, -3, 1, -3, -3},
                        {-inf, -inf, -3, -3, -3, 1, -3},
                        {-inf, -inf, -3, -3, -3, -3, 1}
                    };
                    if (!strcmp(argv[i],"default_gamma"))
                        ++i;
                    else
                        for (int a=0; a<7; ++a)
                            for (int b=0; b<7; ++b)
                                gamma_[a][b]=str2double(argv[i++]);
                    globals.emplace_back(argv[i], gamma_, str2double(argv[i + 1]), str2double(argv[i + 2]), str2double(argv[i + 3]), str2double(argv[i + 4]), str2double(argv[i + 5]), n++, file2browheel[argv[i]]);
                    i += 5;
                    for (auto &node : nodes)
                    {
                        if (node.name == argv[i + 1])
                            globals.back().tail = &node;
                        if (node.name == argv[i + 2])
                            globals.back().head = &node;
                    }
                    i += 2;
                }
                --i;
            }
        std::vector<bool> visit(locals.size() + globals.size());
        RemoveEdge(roots, false, visit, locals, globals);
        RemoveEdge(targets, true, visit, locals, globals);
        SCCTS(locals, globals);
        for (auto &edge : local_crosses)
        {
            edge.tail->out_local_crosses.push_back(&edge);
            edge.head->in_local_crosses.push_back(&edge);
        }
        for (auto &edge : local_circuits)
        {
            edge.tail->out_local_circuits.push_back(&edge);
            edge.head->in_local_circuits.push_back(&edge);
        }
        for (auto &edge : global_crosses)
        {
            edge.tail->out_global_crosses.push_back(&edge);
            edge.head->in_global_crosses.push_back(&edge);
        }
        for (auto &edge : global_circuits)
        {
            edge.tail->out_global_circuits.push_back(&edge);
            edge.head->in_global_circuits.push_back(&edge);
        }
    }

    double str2double(const char *str)
    {
        if (!strcmp(str, "inf"))
            return inf;
        else if (!strcmp(str, "-inf"))
            return -inf;
        else
            return std::stod(str);
    }

    void DeepFirstLine(Node *node, bool reverse, std::vector<bool> &visit, std::list<EdgeLocal> &locals, std::list<EdgeGlobal> &globals)
    {
        for (auto &edge : locals)
            if (reverse ? edge.head == node : edge.tail == node)
                if (!visit[edge.n])
                {
                    visit[edge.n] = true;
                    DeepFirstLine(reverse ? edge.tail : edge.head, reverse, visit, locals, globals);
                }
        for (auto &edge : globals)
            if (reverse ? edge.head == node : edge.tail == node)
                if (!visit[edge.n])
                {
                    visit[edge.n] = true;
                    DeepFirstLine(reverse ? edge.tail : edge.head, reverse, visit, locals, globals);
                }
    }

    void RemoveEdge(std::vector<Node *> &nodes, bool reverse, std::vector<bool> &visit, std::list<EdgeLocal> &locals, std::list<EdgeGlobal> &globals)
    {
        for (auto &edge : locals)
            visit[edge.n] = false;
        for (auto &edge : globals)
            visit[edge.n] = false;
        for (auto &node : nodes)
            DeepFirstLine(node, reverse, visit, locals, globals);
        for (auto iter = locals.begin(); iter != locals.end();)
            if (!visit[iter->n])
                iter=locals.erase(iter);
            else
                ++iter;
        for (auto iter = globals.begin(); iter != globals.end();)
            if (!visit[iter->n])
                iter=globals.erase(iter);
            else
                ++iter;
    }

    void DeepFirst(Node *node, std::stack<Node *> &stack, int &id, std::vector<bool> &visit, std::vector<bool> &in_stack, std::vector<int> &disc, std::vector<int> &low, std::list<EdgeLocal> &locals, std::list<EdgeGlobal> &globals)
    {
        int n = node - nodes.data();
        visit[n] = true;
        disc[n] = id;
        low[n] = id;
        ++id;
        stack.push(node);
        in_stack[n] = true;
        for (auto &edge : locals)
            if (edge.tail == node)
            {
                if (!visit[edge.head - nodes.data()])
                {
                    DeepFirst(edge.head, stack, id, visit, in_stack, disc, low, locals, globals);
                    low[n] = low[n] < low[edge.head - nodes.data()] ? low[n] : low[edge.head - nodes.data()];
                }
                else if (in_stack[edge.head - nodes.data()])
                    low[n] = low[n] < disc[edge.head - nodes.data()] ? low[n] : disc[edge.head - nodes.data()];
            }
        for (auto &edge : globals)
            if (edge.tail == node)
            {
                if (!visit[edge.head - nodes.data()])
                {
                    DeepFirst(edge.head, stack, id, visit, in_stack, disc, low, locals, globals);
                    low[n] = low[n] < low[edge.head - nodes.data()] ? low[n] : low[edge.head - nodes.data()];
                }
                else if (in_stack[edge.head - nodes.data()])
                    low[n] = low[n] < disc[edge.head - nodes.data()] ? low[n] : disc[edge.head - nodes.data()];
            }
        if (disc[n] == low[n])
        {
            sccs.emplace_front();
            Node *top;
            do
            {
                top = stack.top();
                stack.pop();
                in_stack[top - nodes.data()] = false;
                sccs.front().nodes.push_back(top);
            } while (top != node);
            for (auto nnode : sccs.front().nodes)
            {
                nnode->A.resize(sccs.front().nodes.size());
                for (auto &edge : locals)
                    if (edge.tail == nnode)
                    {
                        if (std::find(sccs.front().nodes.begin(), sccs.front().nodes.end(), edge.head) == sccs.front().nodes.end())
                        {
                            local_crosses.emplace_back(edge);
                            sccs.front().local_crosses.push_back(&local_crosses.back());
                        }
                        else
                        {
                            local_circuits.emplace_back(edge);
                            sccs.front().local_circuits.push_back(&local_circuits.back());
                        }
                    }
                for (auto &edge : globals)
                    if (edge.tail == nnode)
                    {
                        if (std::find(sccs.front().nodes.begin(), sccs.front().nodes.end(), edge.head) == sccs.front().nodes.end())
                        {
                            global_crosses.emplace_back(edge);
                            sccs.front().global_crosses.push_back(&global_crosses.back());
                        }
                        else
                        {
                            global_circuits.emplace_back(edge);
                            sccs.front().global_circuits.push_back(&global_circuits.back());
                        }
                    }
            }
        }
    }

    void SCCTS(std::list<EdgeLocal> &locals, std::list<EdgeGlobal> &globals)
    {
        std::vector<bool> visit(nodes.size()), in_stack(nodes.size());
        std::vector<int> disc(nodes.size()), low(nodes.size());
        for (size_t n = 0; n < nodes.size(); ++n)
        {
            visit[n] = false;
            in_stack[n] = false;
        }
        std::stack<Node *> stack;
        int id = 0;
        for (size_t n = 0; n < nodes.size(); ++n)
            if (!visit[n])
                DeepFirst(&nodes[n], stack, id, visit, in_stack, disc, low, locals, globals);
    }

    int draw(std::string file)
    {
        std::ofstream fout(file);
        fout << "digraph align_graph\n{\n";
        fout << "\tnode[shape=circle]\n";
        for (auto root : roots)
            fout << '\t' << root->name << "[shape=square]\n";
        for (auto target : targets)
            fout << '\t' << target->name << "[shape=diamond]\n";
        for (auto &edge : local_crosses)
            fout << '\t' << edge.tail->name << "->" << edge.head->name << "[color=black,label=\"" << edge.name << "\"]\n";
        for (auto &edge : local_circuits)
            fout << '\t' << edge.tail->name << "->" << edge.head->name << "[color=red,label=\"" << edge.name << "\"]\n";
        for (auto &edge : global_crosses)
            fout << '\t' << edge.tail->name << "->" << edge.head->name << "[color=green,label=\"" << edge.name << "\"]\n";
        for (auto &edge : global_circuits)
            fout << '\t' << edge.tail->name << "->" << edge.head->name << "[color=blue,label=\"" << edge.name << "\"]\n";
        fout << "}\n";
        fout.close();
        return system(("dot -Tpdf " + file + " > " + file + ".pdf").c_str());
    }
};

#endif
