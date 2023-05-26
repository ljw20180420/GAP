#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <cfloat>
#include <sstream>
#include "BroWheel.h"
#include <forward_list>
#include <list>

constexpr static const double inf = std::numeric_limits<double>::infinity();

struct Dot // cannot use constructor because on Memory
{
    int n; // n determines the type of Dot
    size_t s;
    int w;
    double val;
    Dot **sources; // apply on memory
    int s_sz;      // s_sz<0 means visited
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

    Edge(std::string name_, std::vector<double> &gamma_, double ve_, double ue_, double vf_, double uf_, double T_, int n_)
    {
        name = name_;
        for (int i = 0, k = 0; i <= sigma; ++i)
            for (int j = 0; j <= sigma; ++j, ++k)
                gamma[i][j] = gamma_[k];
        ve = ve_;
        ue = ue_;
        vf = vf_;
        uf = uf_;
        T = T_;
        n = n_;
    }
};

struct EdgeLocal : Edge
{
    std::string &seq;
    double vfp;
    double ufp;
    double vfm;
    double ufm;

    EdgeLocal(std::string name_, std::vector<double> &gamma_, double ve_, double ue_, double vf_, double uf_, double T_, int n_, std::string &seq_, double vfp_, double ufp_, double vfm_, double ufm_)
        : Edge(name_, gamma_, ve_, ue_, vf_, uf_, T_, n_), seq(seq_)
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
};

struct EdgeGlobal : Edge
{
    BroWheel &browheel;
    std::vector<double> C;
    std::vector<std::list<std::pair<size_t,int>>> Cdelta;
    std::list<TrackNode> tracknodes;

    EdgeGlobal(std::string name_, std::vector<double> &gamma_, double ve_, double ue_, double vf_, double uf_, double T_, int n_, BroWheel &browheel_)
        : Edge(name_, gamma_, ve_, ue_, vf_, uf_, T_, n_), browheel(browheel_)
    {
    }
};

struct EdgeGlobalCross : EdgeGlobal
{
};

struct EdgeGlobalCircuit : EdgeGlobal
{
    std::list<SNC> sncs;
    std::list<std::list<SNC>::iterator> vs;

    Dot *G;
    Dot **D0;
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
    Dot **A;
    Dot *B;
    std::vector<std::map<EdgeGlobalCross*,std::list<std::pair<size_t,int>>>> AdeltaCross;
    std::vector<std::map<EdgeGlobalCircuit*,std::list<std::pair<size_t,int>>>> AdeltaCircuit;

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
    std::vector<EdgeLocalCross> local_crosses;
    std::vector<EdgeLocalCircuit> local_circuits;
    std::vector<EdgeGlobalCross> global_crosses;
    std::vector<EdgeGlobalCircuit> global_circuits;

    std::list<SCC> sccs;

    Graph(int argc, char **argv, std::map<std::string, std::string> &file2seq, std::map<std::string, BroWheel> &file2browheel)
    {
        for (int i = 1; i < argc; ++i)
            if (!strcmp(argv[i], "--nodes"))
            {
                while (argv[++i][0] != '-')
                    nodes.emplace_back(argv[i], std::stod(argv[++i]), std::stod(argv[++i]));
                --i;
            }
        std::vector<EdgeLocal> locals;
        std::vector<EdgeGlobal> globals;
        for (int i = 1, n = 0; i < argc; ++i)
            if (!strcmp(argv[i],"--roots"))
            {
                while (argv[++i][0] != '-')
                    for (auto &node : nodes)
                        if (node.name == argv[i])
                        {
                            roots.push_back(&node);
                            break;
                        }
                --i;
            }
            else if (!strcmp(argv[i],"--targets"))
            {
                while (argv[++i][0] != '-')
                    for (auto &node : nodes)
                        if (node.name == argv[i])
                            targets.push_back(&node);
                --i;
            }
            else if (!strcmp(argv[i],"--locals"))
            {
                while (argv[++i][0] != '-')
                {
                    std::vector<double> gamma_;
                    for (int k = 0; k < (sigma + 1) * (sigma + 1); ++k)
                        gamma_.push_back(std::stod(argv[i++]));
                    locals.emplace_back(argv[i], gamma_, std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), n++, file2seq[argv[++i]], std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]));
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
            else if (!strcmp(argv[i],"--globals"))
            {
                while (argv[++i][0] != '-')
                {
                    std::vector<double> gamma_;
                    for (int k = 0; k < (sigma + 1) * (sigma + 1); ++k)
                        gamma_.push_back(std::stod(argv[i++]));
                    globals.emplace_back(argv[i], gamma_, std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i]), n++, file2browheel[argv[++i]]);
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

    void DeepFirstLine(Node *node, bool reverse, std::vector<bool> &visit, std::vector<EdgeLocal> &locals, std::vector<EdgeGlobal> &globals)
    {
        for (int i = 0; i < locals.size() + globals.size(); ++i)
        {
            Edge *edge = i < locals.size() ? (Edge *)&locals[i] : (Edge *)&globals[i - locals.size()];
            if (reverse ? edge->head == node : edge->tail == node)
                if (!visit[edge->n])
                {
                    visit[edge->n] = true;
                    DeepFirstLine(reverse ? edge->tail : edge->head, reverse, visit, locals, globals);
                }
        }
    }

    void RemoveEdge(std::vector<Node *> &nodes, bool reverse, std::vector<bool> &visit, std::vector<EdgeLocal> &locals, std::vector<EdgeGlobal> &globals)
    {
        for (int i = 0; i < locals.size() + globals.size(); ++i)
        {
            Edge *edge = i < locals.size() ? (Edge *)&locals[i] : (Edge *)&globals[i - locals.size()];
            visit[edge->n] = false;
        }
        for (auto &node : nodes)
            DeepFirstLine(node, reverse, visit, locals, globals);
        for (auto iter = locals.begin(); iter != locals.end(); ++iter)
            if (!visit[iter->n])
                locals.erase(iter);
        for (auto iter = globals.begin(); iter != globals.end(); ++iter)
            if (!visit[iter->n])
                globals.erase(iter);
    }

    void DeepFirst(Node *node, std::stack<Node *> &stack, int &id, std::vector<bool> &visit, std::vector<bool> &in_stack, std::vector<int> &disc, std::vector<int> &low, std::vector<EdgeLocal> &locals, std::vector<EdgeGlobal> &globals)
    {
        int n = node - nodes.data();
        visit[n] = true;
        disc[n] = id;
        low[n] = id;
        ++id;
        stack.push(node);
        in_stack[n] = true;
        for (int i = 0; i < locals.size() + globals.size(); ++i)
        {
            Edge *edge = i < locals.size() ? (Edge *)&locals[i] : (Edge *)&globals[i - locals.size()];
            if (edge->tail == node)
            {
                if (!visit[edge->head - nodes.data()])
                {
                    DeepFirst(edge->head, stack, id, visit, in_stack, disc, low, locals, globals);
                    low[n] = low[n] < low[edge->head - nodes.data()] ? low[n] : low[edge->head - nodes.data()];
                }
                else if (in_stack[edge->head - nodes.data()])
                    low[n] = low[n] < disc[edge->head - nodes.data()] ? low[n] : disc[edge->head - nodes.data()];
            }
        }
        if (disc[n] = low[n])
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

    void SCCTS(std::vector<EdgeLocal> &locals, std::vector<EdgeGlobal> &globals)
    {
        std::vector<bool> visit(nodes.size()), in_stack(nodes.size());
        std::vector<int> disc(nodes.size()), low(nodes.size());
        for (int n=0; n<nodes.size(); ++n)
        {
            visit[n] = false;
            in_stack[n] = false;
        }
        std::stack<Node *> stack;
        int id = 0;
        for (int n=0; n<nodes.size(); ++n)
            if (!visit[n])
                DeepFirst(&nodes[n], stack, id, visit, in_stack, disc, low, locals, globals);
    }

    void draw(std::string file)
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
        system(("dot -Tpdf " + file + " > " + file + ".pdf").c_str());
    }
};

#endif
