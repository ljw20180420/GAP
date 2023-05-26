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

struct Node;

struct Edge
{
    const static int sigma = 6;

    double gamma[sigma + 1][sigma + 1];
    std::pair<std::string,Node *> tail;
    std::pair<std::string,Node *> head;
    double ve;
    double ue;
    double vf;
    double uf;
    double T;
    bool is_global;
    int n;
    SCC* scc_ptr;

    Edge(std::vector<double> gamma_, std::string &tail_name_, std::string &head_name_, double ve_, double ue_, double vf_, double uf_, double T_, bool is_global_, int n_)
    {
        for (int i = 0, k = 0; i <= sigma; ++i)
            for (int j = 0; j <= sigma; ++j, ++k)
                gamma[i][j] = gamma_[k];
        tail.first = tail_name_;
        head.first = head_name_;
        ve = ve_;
        ue = ue_;
        vf = vf_;
        uf = uf_;
        T = T_;
        is_global=is_global_;
        n=n_;
    }
};

struct EdgeLocal : Edge
{
    std::string &seq;
    double vfp;
    double ufp;
    double vfm;
    double ufm;

    EdgeLocal(std::vector<double> gamma_, std::string &tail_name_, std::string &head_name_, double ve_, double ue_, double vf_, double uf_, double T_, bool is_global_, int n_, std::string &seq_, double vfp_, double ufp_, double vfm_, double ufm_)
        : Edge(gamma_, tail_name_, head_name_, ve_, ue_, vf_, uf_, T_, is_global_, n_), seq(seq_)
    {
        vfp = vfp_;
        ufp = ufp_;
        vfm = vfm_;
        ufm = ufm_;
    }
};

template <typename U>
struct EdgeGlobal : Edge
{
    BroWheel<U> &browheel;

    EdgeGlobal(std::vector<double> gamma_, std::string &tail_name_, std::string &head_name_, double ve_, double ue_, double vf_, double uf_, double T_, bool is_global_, int n_, BroWheel<U> &browheel_)
        : Edge(tail_name_, head_name_, ve_, ue_, vf_, uf_, T_, gamma_, is_global_, n_), browheel(browheel_)
    {
    }
};

struct Node
{
    double ve;
    double ue;
    SCC* scc_ptr;
    std::set<std::string> out_edges;
    std::set<std::string> in_edges;

    Node(double ve_, double ue_)
    {
        ve = ve_;
        ue = ue_;
    }
};

struct SCC
{
    std::set<std::string> node_names;
    std::set<std::string> local_names;
    std::set<std::string> global_names;
};

template <typename U>
struct Graph
{
    const static int sigma = 6;

    std::set<std::string> root_names;
    std::set<std::string> target_names;
    std::map<std::string, Node> nodes;
    std::map<std::string, EdgeLocal> locals;
    std::map<std::string, EdgeGlobal<U>> globals;

    std::list<SCC> sccs;

    std::vector<std::string> edge_names;

    Graph(int argc, char **argv, std::map<std::string, std::string> &file2seq, std::map<std::string, BroWheel<U>> &file2browheel)
    {
        for (int i = 1; i < argc; ++i)
            switch (argv[i])
            {
            case "--root_names":
                for (int j = i + 1; j < argc && argv[j][0] != '-'; ++j)
                    root_names.insert(argv[j]);
                break;
            case "--target_names":
                for (int j = i + 1; j < argc && argv[j][0] != '-'; ++j)
                    target_names.insert(argv[j]);
                break;
            case "--nodes":
                for (int j = i + 1; j < argc && argv[j][0] != '-'; ++j)
                    nodes[argv[j]] = Node(std::stod(argv[++j]), std::stod(argv[++j]));
                break;
            case "--locals":
                for (int j = i + 1; j < argc && argv[j][0] != '-'; ++j)
                {
                    std::vector<double> gamma_;
                    for (int k = 0; k < (sigma + 1) * (sigma + 1); ++k)
                        gamma_.push_back(std::stod(argv[j++]));
                    edge_names.push_back(argv[j]);
                    locals[argv[j]] = EdgeLocal(gamma_, argv[++j], argv[++j], std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), false, edge_names.size()-1, file2seq[argv[++j]], std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]));
                }
                break;
            case "--globals":
                for (int j = i + 1; j < argc && argv[j][0] != '-'; ++j)
                {
                    std::vector<double> gamma_;
                    for (int k = 0; k < (sigma + 1) * (sigma + 1); ++k)
                        gamma_.push_back(std::stod(argv[j++]));
                    edge_names.push_back(argv[j]);
                    globals[argv[j]] = EdgeGlobal<U>(gamma_, argv[++j], argv[++j], std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), std::stod(argv[++j]), true, edge_names.size()-1, file2browheel[argv[++j]]);
                }
                break;
            }
        for (auto &local : locals)
        {
            nodes[local.second.tail.first].out_edges.insert(local.first);
            nodes[local.second.head.first].in_edges.insert(local.first);
        }
        for (auto &global : globals)
        {
            nodes[global.second.tail.first].out_edges.insert(global.first);
            nodes[global.second.head.first].in_edges.insert(global.first);
        }
        RemoveEdge(root_names, false);
        RemoveEdge(target_names, true);
        SCCTS();
        for (auto &scc : sccs)
        {
            for (auto &node_name : scc.node_names)
                nodes[node_name].scc_ptr=&scc;
            for (auto &local_name : scc.local_names)
                locals[local_name].scc_ptr=&scc;
            for (auto &global_name : scc.global_names)
                globals[global_name].scc_ptr=&scc;
        }
    }

    void DeepFirstLine(std::string &node_name, bool reverse, std::map<std::string, bool> &visit)
    {
        for (auto &edge_name : (reverse ? nodes[node_name].in_edges : nodes[node_name].out_edges))
        {
            auto &edge = globals.find(edge_name)==globals.end() ? locals[edge_name] : globals[edge_name];
            if (!edge.visit)
            {
                edge.visit = true;
                DeepFirstLine(reverse ? edge.tail.first : edge.head.first, reverse, visit);
            }
        }
    }

    void RemoveEdge(std::set<std::string> &node_names, bool reverse)
    {
        std::map<std::string, bool> visit;
        for (auto &edge_name: edge_names)
            visit[edge_name] = false;
        for (auto &node_name : node_names)
            DeepFirstLine(node_name, reverse, visit);
        for (auto &local : locals)
            if (!visit[local.first])
                locals.erase(local.first);
        for (auto &global : globals)
            if (!visit[global.first])
                globals.erase(global.first);
    }

    void DeepFirst(std::string &node_name, std::stack<std::string> &stack_names, int &id, std::map<std::string, bool> &visit, std::map<std::string, bool> &in_stack, std::map<std::string, int> &disc, std::map<std::string, int> &low)
    {
        visit[node_name] = true;
        disc[node_name] = id;
        low[node_name] = id;
        ++id;
        stack_names.push(node_name);
        in_stack[node_name] = true;
        for (auto &edge_name : nodes[node_name].out_edges)
        {
            auto &edge = globals.find(edge_name)==globals.end() ? locals[edge_name] : globals[edge_name];
            if (!visit[edge.head.first])
            {
                DeepFirst(edge.head.first, stack_names, id, visit, in_stack, disc, low);
                low[node_name] = low[node_name] < low[edge.head.first] ? low[node_name] : low[edge.head.first];
            }
            else if (in_stack[edge.head.first])
                low[node_name] = low[node_name] < disc[edge.head.first] ? low[node_name] : disc[edge.head.first];
        }
        if (disc[node_name] = low[node_name])
        {
            sccs.emplace_front();
            std::string top_name;
            do
            {
                top_name = stack_names.top();
                stack_names.pop();
                in_stack[top_name] = false;
                sccs.front().node_names.insert(top_name);
            } while (top_name != node_name);
            sccs.emplace(std::next(sccs.begin()));
            for (auto &node_name : sccs.front().node_names)
                for (auto &edge_name : nodes[node_name].out_edges)
                {
                    auto &edge = globals.find(edge_name)==globals.end() ? locals[edge_name] : globals[edge_name];
                    auto it = sccs.begin();
                    if (it->node_names.find(edge.head.first) == it->node_names.end())
                        ++it;
                    auto &edge_names = globals.find(edge_name)==globals.end() ? it->local_names : it->global_names;
                    edge_names.insert(edge_name);
                }
        }
    }

    void SCCTS()
    {
        std::map<std::string, bool> visit, in_stack;
        std::map<std::string, int> disc, low;
        for (auto &node : nodes)
        {
            visit[node.first] = false;
            in_stack[node.first] = false;
        }
        std::stack<std::string> stack_names;
        int id = 0;
        for (auto &node : nodes)
            if (!visit[node.first])
                DeepFirst(node.first, stack_names, id, visit, in_stack, disc, low);
    }

    void draw(std::string file)
    {
        std::ofstream fout(file);
        fout << "digraph align_graph\n{\n";
        fout << "\tnode[shape=circle]\n";
        for (auto &root_name : root_names)
            fout << '\t' << root_name << "[shape=square]\n";
        for (auto &target_name : target_names)
            fout << '\t' << target_name << "[shape=diamond]\n";
        for (auto &local : locals)
            fout << '\t' << local.second.tail.first << "->" << local.second.head.first << "[color=black,label=\"" << local.first << "\"]\n";
        for (auto &global : globals)
            fout << '\t' << global.second.tail.first << "->" << global.second.head.first << "[color=red,label=\"" << global.first << "\"]\n";
        fout << "}\n";
        fout.close();
        system(("dot -Tpdf " + file + " > " + file + ".pdf").c_str());
    }
};

#endif
