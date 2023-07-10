#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <cfloat>
#include <list>
#include "BroWheel.h"
#include <boost/program_options.hpp>

constexpr static const double inf = std::numeric_limits<double>::infinity();

struct Dot
{
    static const int DotQ = -3, DotAbar = -2, DotB = -1;

    int n; // n determines the type of Dot
    int64_t s;
    int w;
    double val;
    int64_t fs; // first source
    int s_sz;
    int lambda;
    int64_t id;

    static int nidx_trans(int nidx)
    {
        return -nidx - 4;
    }
};

struct Node;

struct Edge
{
    const static int sigma = 6;

    int n;
    std::string name;
    double gamma[sigma + 1][sigma + 1] = {
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -3, -3, -3, -3, -3},
        {-inf, -inf, -3, 1, -3, -3, -3},
        {-inf, -inf, -3, -3, 1, -3, -3},
        {-inf, -inf, -3, -3, -3, 1, -3},
        {-inf, -inf, -3, -3, -3, -3, 1}};
    double ve, ue, vf, uf, T, min_score;
    Node *tail;
    Node *head;
};

struct EdgeLocal : Edge
{
    double vfp;
    double ufp;
    double vfm;
    double ufm;
    std::vector<double> gf;
    double gfm;
    std::vector<double> gfpT;
};

template <typename T>
void type_initial(T *&Ts, std::initializer_list<T **> pTss, std::initializer_list<std::unique_ptr<T *[]> *> pTsss, int *Ss, uint64_t Omax)
{
    for (T **pTs : pTss)
    {
        *pTs = Ts;
        Ts += Omax + 1;
    }
    int si = 0;
    for (auto it = pTsss.begin(); it != pTsss.end(); ++it, ++si)
    {
        (*it)->reset(new T *[Ss[si] + 1]);
        for (int s = 0; s <= Ss[si]; ++s, Ts += Omax + 1)
            (**it)[s] = Ts;
    }
}

void extend_ss_ns(int pTss_sz, int n, int *&fps, int *&fpn)
{
    for (int i = 0; i < pTss_sz; ++i, ++fps, ++fpn)
    {
        *fps = 0;
        *fpn = n;
    }
}

void extend_ss_ns(int pTsss_sz, int *Ss, int n, int *&fps, int *&fpn)
{
    for (int i = 0; i < pTsss_sz; ++i)
        for (int j = 0; j <= Ss[i]; ++j, ++fps, ++fpn)
        {
            *fps = j;
            *fpn = n;
        }
}

struct EdgeLocalCross;
struct EdgeLocalCircuit;
struct EdgeGlobalCircuit;

struct Node
{
    bool is_root, is_target;
    int n, scc_id, scc_sz;
    std::string name;
    double ve, ue;

    std::unique_ptr<double *[]> Avals;
    double *Bvals, *pAbarval;
    std::unique_ptr<double ***[]> Asourcess;
    double ***Bsourcess;
    std::unique_ptr<int *[]> As_szs;
    int *Bs_szs;
    std::unique_ptr<int64_t *[]> Aids;
    int64_t *Bids, *pAbarid;

    std::unique_ptr<std::deque<double *>[]> AdeltaDot;

    struct GlobalSuffix
    {
        Edge *edge;
        int64_t start;
        int lambda;

        GlobalSuffix(Edge *edge_, int64_t start_, int lambda_)
            : edge(edge_), start(start_), lambda(lambda_)
        {
        }
    };

    std::unique_ptr<std::deque<GlobalSuffix>[]> AdeltaGlobal;

    int64_t get_tsn(uint64_t Aso, uint64_t Omax)
    {
        return (Omax + 1) * (2 + Aso * (scc_sz - 1));
    }

    int *get_Ss()
    {
        return new int[1]{scc_sz - 1};
    }

    int64_t *get_SE(int Omax)
    {
        return new int64_t[4]{0, Omax + 1, 2 * (Omax + 1), int64_t(Omax + 1) * (1 + scc_sz)};
    }

    int *get_steps(uint64_t Aso)
    {
        return new int[3]{2, 0, Aso};
    }

    void apply_memory(uint64_t Aso, uint64_t Omax, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn, double *&rpval, int64_t *&rpid)
    {
        std::unique_ptr<int[]> Ss(new int[1]{scc_sz - 1});
        extend_ss_ns(1, -1, fps, fpn);
        extend_ss_ns(1, Ss.get(), Dot::nidx_trans(n), fps, fpn);
        type_initial(fpval, {&Bvals}, {&Avals}, Ss.get(), Omax);
        pAbarval = --rpval;
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {&Bsourcess}, {&Asourcess}, Ss.get(), Omax);
        std::unique_ptr<int64_t[]> SE(get_SE(Omax));
        std::unique_ptr<int[]> steps(get_steps(Aso));
        for (int i = 0; i < 3; ++i)
            for (int64_t j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {&Bs_szs}, {&As_szs}, Ss.get(), Omax);
        type_initial(fpid, {&Bids}, {&Aids}, Ss.get(), Omax);
        pAbarid = --rpid;

        AdeltaDot.reset(new std::deque<double *>[Omax + 1]);
        AdeltaGlobal.reset(new std::deque<Node::GlobalSuffix>[Omax + 1]);
    }

    void updateA0(int w, double source_val, double *source)
    {
        if (source_val >= Avals[0][w])
        {
            if (source_val > Avals[0][w])
            {
                Avals[0][w] = source_val;
                AdeltaDot[w].clear();
                AdeltaGlobal[w].clear();
            }
            AdeltaDot[w].emplace_back(source);
        }
    }

    void updateA0(int w, double source_val, Edge *edge_, int64_t start_, int lambda_)
    {
        if (source_val >= Avals[0][w])
        {
            if (source_val > Avals[0][w])
            {
                Avals[0][w] = source_val;
                AdeltaDot[w].clear();
                AdeltaGlobal[w].clear();
            }
            AdeltaGlobal[w].emplace_back(edge_, start_, lambda_);
        }
    }

    void clearAdelta(int W)
    {
        for (int w = 0; w <= W; ++w)
        {
            AdeltaDot[w].clear();
            AdeltaGlobal[w].clear();
        }
    }
};

struct EdgeLocalCross : EdgeLocal
{
    std::unique_ptr<double *[]> Evals, Fvals, Gvals;
    std::unique_ptr<double ***[]> Esourcess, Fsourcess, Gsourcess;
    std::unique_ptr<int *[]> Es_szs, Fs_szs, Gs_szs;
    std::unique_ptr<int64_t *[]> Eids, Fids, Gids;

    EdgeLocalCross(EdgeLocal &edge)
        : EdgeLocal(edge)
    {
    }

    int get_trn(uint64_t ref_sz)
    {
        return 3 * (ref_sz + 1);
    }

    int64_t get_tnn(uint64_t Omax, uint64_t ref_sz)
    {
        return (Omax + 1) * get_trn(ref_sz);
    }

    int64_t get_tsn(uint64_t Omax, uint64_t ref_sz)
    {
        return (Omax + 1) * 8 * (ref_sz + 1);
    }

    int *get_Ss(uint64_t ref_sz)
    {
        int S = ref_sz;
        return new int[3]{S, S, S};
    }

    int64_t *get_SE(uint64_t Omax, uint64_t ref_sz)
    {
        int64_t mid = (Omax + 1) * 2 * (ref_sz + 1);
        return new int64_t[3]{0, mid, get_tnn(Omax, ref_sz)};
    }

    int *get_steps()
    {
        return new int[2]{2, 4};
    }

    void apply_memory(uint64_t Omax, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn, uint64_t ref_sz)
    {
        std::unique_ptr<int[]> Ss(get_Ss(ref_sz));
        extend_ss_ns(3, Ss.get(), n, fps, fpn);
        type_initial(fpval, {}, {&Evals, &Fvals, &Gvals}, Ss.get(), Omax);
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {}, {&Esourcess, &Fsourcess, &Gsourcess}, Ss.get(), Omax);
        std::unique_ptr<int64_t[]> SE(get_SE(Omax, ref_sz));
        std::unique_ptr<int[]> steps(get_steps());
        for (int i = 0; i < 2; ++i)
            for (int64_t j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {}, {&Es_szs, &Fs_szs, &Gs_szs}, Ss.get(), Omax);
        type_initial(fpid, {}, {&Eids, &Fids, &Gids}, Ss.get(), Omax);
    }
};

struct EdgeLocalCircuit : EdgeLocal
{
    std::unique_ptr<double *[]> Evals, F0vals, G0vals, Gvals, D0vals, DXvals;
    double *Dvals;
    std::unique_ptr<double ***[]> Esourcess, F0sourcess, G0sourcess, Gsourcess;
    double ***Dsourcess;
    std::unique_ptr<int *[]> Es_szs, F0s_szs, G0s_szs, Gs_szs;
    int *Ds_szs;
    std::unique_ptr<int64_t *[]> Eids, F0ids, G0ids, Gids, D0ids, DXids;
    int64_t *Dids;
    std::unique_ptr<uint8_t []> DXbits;

    EdgeLocalCircuit(EdgeLocal &edge)
        : EdgeLocal(edge)
    {
    }

    int get_trn(uint64_t ref_sz)
    {
        return 1 + 4 * (ref_sz + 1) + 2 * (tail->scc_sz - 1);
    }

    int64_t get_tnn(uint64_t Omax, uint64_t ref_sz)
    {
        return (Omax + 1) * get_trn(ref_sz);
    }

    int64_t get_tsn(uint64_t Omax, uint64_t ref_sz)
    {
        return (Omax + 1) * (1 + 10 * (ref_sz + 1) + 3 * (tail->scc_sz - 1));
    }

    int *get_Ss(uint64_t ref_sz)
    {
        int S = ref_sz;
        int scc_sz = tail->scc_sz;
        return new int[6]{S, S, S, S, scc_sz - 2, scc_sz - 2};
    }

    int64_t *get_SE(uint64_t Omax, uint64_t ref_sz)
    {
        int S = ref_sz;
        int scc_sz = tail->scc_sz;
        return new int64_t[6]{0, Omax + 1, (Omax + 1) * (1 + 2 * (S + 1)), (Omax + 1) * (1 + 4 * (S + 1))};
    }

    int *get_steps()
    {
        return new int[5]{1, 2, 3, 1, 2};
    }

    void apply_memory(uint64_t Omax, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn, uint64_t ref_sz)
    {
        DXbits.reset(new uint8_t[(tail->scc_sz-1)*(Omax+1)]);
        std::unique_ptr<int[]> Ss(get_Ss(ref_sz));
        extend_ss_ns(1, n, fps, fpn);
        extend_ss_ns(6, Ss.get(), n, fps, fpn);
        type_initial(fpval, {&Dvals}, {&Evals, &F0vals, &G0vals, &Gvals, &D0vals, &DXvals}, Ss.get(), Omax);
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {&Dsourcess}, {&Esourcess, &F0sourcess, &G0sourcess, &Gsourcess}, Ss.get(), Omax);
        fpsources += (tail->scc_sz - 1) * 2 * (Omax + 1);
        std::unique_ptr<int64_t[]> SE(get_SE(Omax, ref_sz));
        std::unique_ptr<int[]> steps(get_steps());
        for (int i = 0; i < 3; ++i)
            for (int64_t j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {&Ds_szs}, {&Es_szs, &F0s_szs, &G0s_szs, &Gs_szs}, Ss.get(), Omax);
        fps_sz += (tail->scc_sz - 1) * 2 * (Omax + 1);
        type_initial(fpid, {&Dids}, {&Eids, &F0ids, &G0ids, &Gids, &D0ids, &DXids}, Ss.get(), Omax);
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
    int64_t sr1;
    int64_t sr2;
    int64_t cid; // 0 means not try to add child yet
    SNC *itp;
    SNC *jump;

    SNC(double E_, double F0_, double G0_, double hatG_, int8_t c_, int8_t idl_, int lambda_, int64_t sr1_, int64_t sr2_, int64_t cid_, SNC *itp_, SNC *jump_)
        : E(E_), F0(F0_), G0(G0_), hatG(hatG_), c(c_), idl(idl_), lambda(lambda_), sr1(sr1_), sr2(sr2_), cid(cid_), itp(itp_), jump(jump_)
    {
    }
};

struct EdgeGlobalCircuit : Edge
{
    std::deque<SNC> sncs;
    std::unique_ptr<double *[]> D0vals;
    std::unique_ptr<int64_t *[]> D0ids;

    EdgeGlobalCircuit(Edge &edge)
        : Edge(edge)
    {
    }

    int get_trn()
    {
        return tail->scc_sz - 1;
    }

    int64_t get_tnn(uint64_t Omax)
    {
        return (Omax + 1) * get_trn();
    }

    int *get_Ss()
    {
        return new int[1]{tail->scc_sz - 2};
    }

    void apply_memory(uint64_t Omax, double *&fpval, int64_t *&fpid, int *&fps, int *&fpn)
    {
        std::unique_ptr<int[]> Ss(get_Ss());
        extend_ss_ns(1, Ss.get(), n, fps, fpn);
        type_initial(fpval, {}, {&D0vals}, Ss.get(), Omax);
        type_initial(fpid, {}, {&D0ids}, Ss.get(), Omax);
    }
};

struct SCC
{
    std::vector<Node *> nodes;
    std::vector<EdgeLocalCross *> local_crosses;
    std::vector<EdgeLocalCircuit *> local_circuits;
    std::vector<Edge *> global_crosses;
    std::vector<EdgeGlobalCircuit *> global_circuits;
};

struct Graph
{
    const static int sigma = 6;

    std::deque<Node> nodes;
    std::deque<EdgeLocalCross> local_crosses; // cannot be changed to vector, because emplace_back of vector may change the addresses of elements
    std::deque<EdgeLocalCircuit> local_circuits;
    std::deque<Edge> global_crosses;
    std::deque<EdgeGlobalCircuit> global_circuits;
    std::deque<Edge *> edges;

    std::deque<SCC> sccs;

    std::unique_ptr<double[]> vals;
    std::unique_ptr<double *[]> sources;
    std::unique_ptr<double **[]> sourcess;
    std::unique_ptr<int[]> s_szs;
    std::unique_ptr<int64_t[]> ids;
    std::unique_ptr<int[]> ss, ns;
    int trn;
    int64_t tnn, tsn;

    double *pQval;
    double ***pQsources;
    int *pQs_sz;
    int64_t *pQid;

    boost::program_options::variables_map &vm;
    std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2short;
    std::map<std::string, RankVec> &file2rankvec;
    std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2long;
    std::map<std::string, std::ifstream> file2SA;

    void get_affine(std::vector<double> &g, int seqlen, double initial, double u, double v)
    {
        g.resize(seqlen + 1);
        g[0] = initial;
        for (int s = 0; s < seqlen; ++s)
            if (s == 0)
                g[s + 1] = g[s] + v;
            else
                g[s + 1] = g[s] + u;
    }

    Graph(boost::program_options::variables_map &vm_, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2short_, std::map<std::string, RankVec> &file2rankvec_, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2long_) : vm(vm_), file2short(file2short_), file2rankvec(file2rankvec_), file2long(file2long_)
    {
        for (std::string file : vm["longs"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2SA.count(file))
                continue;
            file2SA[file].open(file+".sa");
        }

        for (std::string nodeinfo : vm["nodes"].as<std::vector<std::string>>())
        {
            std::stringstream ss(nodeinfo);
            std::string name, is_root, is_target, v, u;
            std::getline(std::getline(std::getline(std::getline(std::getline(ss, name, ','), is_root, ','), is_target, ','), v, ','), u, ',');
            
            Node node;
            node.n = nodes.size();
            node.name = name;
            node.is_root = std::stoi(is_root);
            node.is_target = std::stoi(is_target);
            try{node.ve = std::stod(v);}catch(...){node.ve = 0.0;}
            try{node.ue = std::stod(u);}catch(...){node.ue = 0.0;}
            nodes.push_back(std::move(node));
        }

        std::list<EdgeLocal> locals;
        if (vm.count("shorts"))
            for (std::string shortinfo : vm["shorts"].as<std::vector<std::string>>())
            {
                std::stringstream ss(shortinfo);
                std::string name, tail, head, match, mismatch, v, u, vb, ub, T, min_score;
                std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(ss, name, ','), tail, ','), head, ','), match, ','), mismatch, ','), v, ','), u, ','), T, ','), min_score, ','), vb, ','), ub, ',');

                EdgeLocal local;
                local.n = locals.size();
                for (int a = 2; a < 7; ++a)
                    for (int b = 2; b < 7; ++b)
                        if (a==b && a>2)
                            try{local.gamma[a][b]=std::stod(match);}catch(...){local.gamma[a][b]=1.0;}
                        else
                            try{local.gamma[a][b]=std::stod(mismatch);}catch(...){local.gamma[a][b]=-3.0;}
                local.name = name;
                for (Node &node : nodes)
                {
                    if (node.name == head)
                        local.head = &node;
                    if (node.name == tail)
                        local.tail = &node;
                }
                try{local.ve = std::stod(v);}catch(...){local.ve = -5.0;}
                local.vf = local.ve;
                try{local.ue = std::stod(u);}catch(...){local.ue = -2.0;}
                local.uf = local.ue;
                try{local.T = std::stod(T);}catch(...){local.T = -10.0;}
                try{local.min_score = std::stod(min_score);}catch(...){local.min_score = 20.0;}
                try{local.vfp = std::stod(vb);}catch(...){local.vfp = 0.0;}
                local.vfm = local.vfp;
                try{local.ufp = std::stod(ub);}catch(...){local.ufp = 0.0;}
                local.ufm = local.ufp;
                uint64_t ref_sz = file2short[local.name].second;
                get_affine(local.gf, ref_sz, 0, local.uf, local.vf);
                local.gfm = local.vfm + (ref_sz - 1) * local.ufm;
                get_affine(local.gfpT, ref_sz, local.T, local.ufp, local.vfp);
                locals.push_back(local);
            }


        std::list<Edge> globals;
        if (vm.count("longs"))
            for (std::string longinfo : vm["longs"].as<std::vector<std::string>>())
            {
                std::stringstream ss(longinfo);
                std::string name, tail, head, match, mismatch, v, u, vb, ub, T, min_score;
                std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(std::getline(ss, name, ','), tail, ','), head, ','), match, ','), mismatch, ','), v, ','), u, ','), T, ','), min_score, ',');

                Edge global;
                global.n = locals.size() + globals.size();
                for (int a = 2; a < 7; ++a)
                    for (int b = 2; b < 7; ++b)
                        if (a==b && a>2)
                            try{global.gamma[a][b]=std::stod(match);}catch(...){global.gamma[a][b]=1.0;}
                        else
                            try{global.gamma[a][b]=std::stod(mismatch);}catch(...){global.gamma[a][b]=-3.0;}
                global.name = name;
                for (Node &node : nodes)
                {
                    if (node.name == head)
                        global.head = &node;
                    if (node.name == tail)
                        global.tail = &node;
                }
                try{global.ve = std::stod(v);}catch(...){global.ve = -5.0;}
                global.vf = global.ve;
                try{global.ue = std::stod(u);}catch(...){global.ue = -2.0;}
                global.uf = global.ue;
                try{global.T = std::stod(T);}catch(...){global.T = -10.0;}
                try{global.min_score = std::stod(min_score);}catch(...){global.min_score = 20.0;}

                globals.push_back(global);
            }

        std::vector<bool> visit(locals.size() + globals.size());
        RemoveEdge(nodes, false, visit, locals, globals);
        RemoveEdge(nodes, true, visit, locals, globals);
        SCCTS(locals, globals);
        for (uint64_t i = 0; i < sccs.size(); ++i)
            for (auto node : sccs[i].nodes)
            {
                node->scc_id = i;
                node->scc_sz = sccs[i].nodes.size();
            }
        edges.resize(locals.size() + globals.size());
        for (auto &edge : locals)
            if (edge.tail->scc_id != edge.head->scc_id)
            {
                local_crosses.emplace_back(edge);
                sccs[edge.tail->scc_id].local_crosses.emplace_back(&local_crosses.back());
                edges[edge.n] = &local_crosses.back();
            }
            else
            {
                local_circuits.emplace_back(edge);
                sccs[edge.tail->scc_id].local_circuits.emplace_back(&local_circuits.back());
                edges[edge.n] = &local_circuits.back();
            }
        for (auto &edge : globals)
            if (edge.tail->scc_id != edge.head->scc_id)
            {
                global_crosses.emplace_back(edge);
                sccs[edge.tail->scc_id].global_crosses.emplace_back(&global_crosses.back());
                edges[edge.n] = &global_crosses.back();
            }
            else
            {
                global_circuits.emplace_back(edge);
                sccs[edge.tail->scc_id].global_circuits.emplace_back(&global_circuits.back());
                edges[edge.n] = &global_circuits.back();
            }
    }

    double str2double(std::string &str)
    {
        return str2double(str.c_str());
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

    void DeepFirstLine(Node *node, bool reverse, std::vector<bool> &visit, std::list<EdgeLocal> &locals, std::list<Edge> &globals)
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

    void RemoveEdge(std::deque<Node> &nodes, bool reverse, std::vector<bool> &visit, std::list<EdgeLocal> &locals, std::list<Edge> &globals)
    {
        for (auto &edge : locals)
            visit[edge.n] = false;
        for (auto &edge : globals)
            visit[edge.n] = false;
        for (auto &node : nodes)
            if ((!reverse && node.is_root) || (reverse && node.is_target))
                DeepFirstLine(&node, reverse, visit, locals, globals);
        for (auto iter = locals.begin(); iter != locals.end();)
            if (!visit[iter->n])
                iter = locals.erase(iter);
            else
                ++iter;
        for (auto iter = globals.begin(); iter != globals.end();)
            if (!visit[iter->n])
                iter = globals.erase(iter);
            else
                ++iter;
    }

    void DeepFirst(Node *node, std::stack<Node *> &stack, int &id, std::vector<bool> &visit, std::vector<bool> &in_stack, std::vector<int> &disc, std::vector<int> &low, std::list<EdgeLocal> &locals, std::list<Edge> &globals)
    {
        int n = node->n;
        visit[n] = true;
        disc[n] = id;
        low[n] = id;
        ++id;
        stack.push(node);
        in_stack[n] = true;
        for (auto &edge : locals)
            if (edge.tail == node)
            {
                if (!visit[edge.head->n])
                {
                    DeepFirst(edge.head, stack, id, visit, in_stack, disc, low, locals, globals);
                    low[n] = low[n] < low[edge.head->n] ? low[n] : low[edge.head->n];
                }
                else if (in_stack[edge.head->n])
                    low[n] = low[n] < disc[edge.head->n] ? low[n] : disc[edge.head->n];
            }
        for (auto &edge : globals)
            if (edge.tail == node)
            {
                if (!visit[edge.head->n])
                {
                    DeepFirst(edge.head, stack, id, visit, in_stack, disc, low, locals, globals);
                    low[n] = low[n] < low[edge.head->n] ? low[n] : low[edge.head->n];
                }
                else if (in_stack[edge.head->n])
                    low[n] = low[n] < disc[edge.head->n] ? low[n] : disc[edge.head->n];
            }
        if (disc[n] == low[n])
        {
            sccs.emplace_front();
            Node *top;
            do
            {
                top = stack.top();
                stack.pop();
                in_stack[top->n] = false;
                sccs.front().nodes.push_back(top);
            } while (top != node);
        }
    }

    void SCCTS(std::list<EdgeLocal> &locals, std::list<Edge> &globals)
    {
        std::vector<bool> visit(nodes.size()), in_stack(nodes.size());
        std::vector<int> disc(nodes.size()), low(nodes.size());
        for (uint64_t n = 0; n < nodes.size(); ++n)
        {
            visit[n] = false;
            in_stack[n] = false;
        }
        std::stack<Node *> stack;
        int id = 0;
        for (uint64_t n = 0; n < nodes.size(); ++n)
            if (!visit[n])
                DeepFirst(&nodes[n], stack, id, visit, in_stack, disc, low, locals, globals);
    }

    int draw(std::string file)
    {
        std::ofstream fout(file);
        fout << "digraph align_graph\n{\n";
        fout << "\tnode[shape=circle]\n";
        for (Node &node : nodes)
            if (node.is_root)
                fout << '\t' << node.name << "[shape=square]\n";
            else if (node.is_target)
                fout << '\t' << node.name << "[shape=diamond]\n";            
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
