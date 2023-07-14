#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <cfloat>
#include <list>
#include "BroWheel.h"
#include <boost/program_options.hpp>

using DOTTYPE = int16_t;
using QUERYSIZE = uint32_t;
using SHORTSIZE = uint32_t;
using IDTYPE = uint64_t;
// using SCORETYPE = int32_t;
using SCORETYPE = double;
using SOURCESIZE = uint32_t;
constexpr static const SCORETYPE inf = std::numeric_limits<SCORETYPE>::infinity();

struct Dot
{
    static const DOTTYPE DotQ = -3, DotAbar = -2, DotB = -1;

    DOTTYPE n; // n determines the type of Dot
    SIZETYPE s;
    QUERYSIZE w;
    SCORETYPE val;
    SIZETYPE fs; // first source
    SOURCESIZE s_sz;
    QUERYSIZE lambda;
    IDTYPE id;

    static DOTTYPE nidx_trans(DOTTYPE nidx)
    {
        return -nidx - 4;
    }
};

struct Node;

struct Edge
{
    DOTTYPE n;
    std::string name;
    SCORETYPE gamma[RankVec::sigma][RankVec::sigma] = {
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -3, -3, -3, -3, -3},
        {-inf, -inf, -3, 1, -3, -3, -3},
        {-inf, -inf, -3, -3, 1, -3, -3},
        {-inf, -inf, -3, -3, -3, 1, -3},
        {-inf, -inf, -3, -3, -3, -3, 1}};
    SCORETYPE ve, ue, vf, uf, T, min_score;
    Node *tail;
    Node *head;
};

struct EdgeLocal : Edge
{
    SCORETYPE vfp, ufp, vfm, ufm, gfm;
    std::vector<SCORETYPE> gf, gfpT;
};

template <typename T>
void type_initial(T *&Ts, std::initializer_list<T **> pTss, std::initializer_list<std::unique_ptr<T *[]> *> pTsss, int64_t *Ss, QUERYSIZE Omax)
{
    for (T **pTs : pTss)
    {
        *pTs = Ts;
        Ts += Omax + 1;
    }
    SIZETYPE si = 0;
    for (typename std::initializer_list<std::unique_ptr<T *[]> *>::iterator it = pTsss.begin(); it != pTsss.end(); ++it, ++si)
    {
        (*it)->reset(new T *[Ss[si] + 1]);
        for (int64_t s = 0; s <= Ss[si]; ++s, Ts += Omax + 1)
            (**it)[s] = Ts;
    }
}

void extend_ss_ns(SIZETYPE pTss_sz, DOTTYPE n, SHORTSIZE *&fps, DOTTYPE *&fpn)
{
    for (SIZETYPE i = 0; i < pTss_sz; ++i, ++fps, ++fpn)
    {
        *fps = 0;
        *fpn = n;
    }
}

void extend_ss_ns(SIZETYPE pTsss_sz, int64_t *Ss, DOTTYPE n, SHORTSIZE *&fps, DOTTYPE *&fpn)
{
    for (SIZETYPE i = 0; i < pTsss_sz; ++i)
        for (int64_t j = 0; j <= Ss[i]; ++j, ++fps, ++fpn)
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
    DOTTYPE n;
    SIZETYPE scc_id, scc_sz;
    std::string name;
    SCORETYPE ve, ue;

    std::unique_ptr<SCORETYPE *[]> Avals;
    SCORETYPE *Bvals, *pAbarval;
    std::unique_ptr<SCORETYPE ***[]> Asourcess;
    SCORETYPE ***Bsourcess;
    std::unique_ptr<SOURCESIZE *[]> As_szs;
    SOURCESIZE *Bs_szs;
    std::unique_ptr<IDTYPE *[]> Aids;
    IDTYPE *Bids, *pAbarid;

    std::unique_ptr<std::deque<SCORETYPE *> []> AdeltaDot;

    struct GlobalSuffix
    {
        Edge *edge;
        SIZETYPE start;
        QUERYSIZE lambda;

        GlobalSuffix(Edge *edge_, SIZETYPE start_, QUERYSIZE lambda_)
            : edge(edge_), start(start_), lambda(lambda_)
        {
        }
    };

    std::unique_ptr<std::deque<GlobalSuffix> []> AdeltaGlobal;

    SIZETYPE get_tsn(SIZETYPE Aso, QUERYSIZE Omax)
    {
        return (2 + Aso * (scc_sz - 1)) * (Omax + 1);
    }

    SIZETYPE *get_SE(QUERYSIZE Omax)
    {
        return new SIZETYPE[4]{0, Omax + 1, 2 * (Omax + 1), (1 + scc_sz) * (Omax + 1)};
    }

    SIZETYPE *get_steps(SIZETYPE Aso)
    {
        return new SIZETYPE[3]{2, 0, Aso};
    }

    void apply_memory(SIZETYPE Aso, QUERYSIZE Omax, SCORETYPE *&fpval, SCORETYPE **&fpsource, SCORETYPE ***&fpsources, SOURCESIZE *&fps_sz, IDTYPE *&fpid, SHORTSIZE *&fps, DOTTYPE *&fpn, SCORETYPE *&rpval, IDTYPE *&rpid)
    {
        std::unique_ptr<int64_t []> Ss(new int64_t[1]{scc_sz - 1});
        extend_ss_ns(1, -1, fps, fpn);
        extend_ss_ns(1, Ss.get(), Dot::nidx_trans(n), fps, fpn);
        type_initial(fpval, {&Bvals}, {&Avals}, Ss.get(), Omax);
        pAbarval = --rpval;
        SCORETYPE ***fpsources_old = fpsources;
        type_initial(fpsources, {&Bsourcess}, {&Asourcess}, Ss.get(), Omax);
        std::unique_ptr<SIZETYPE []> SE(get_SE(Omax));
        std::unique_ptr<SIZETYPE []> steps(get_steps(Aso));
        for (SIZETYPE i = 0; i < 3; ++i)
            for (SIZETYPE j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {&Bs_szs}, {&As_szs}, Ss.get(), Omax);
        type_initial(fpid, {&Bids}, {&Aids}, Ss.get(), Omax);
        pAbarid = --rpid;

        AdeltaDot.reset(new std::deque<SCORETYPE *>[Omax + 1]);
        AdeltaGlobal.reset(new std::deque<Node::GlobalSuffix>[Omax + 1]);
    }

    void updateA0(QUERYSIZE w, SCORETYPE source_val, SCORETYPE *source)
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

    void updateA0(QUERYSIZE w, SCORETYPE source_val, Edge *edge_, SIZETYPE start_, QUERYSIZE lambda_)
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

    void clearAdelta(QUERYSIZE W)
    {
        for (SIZETYPE w = 0; w <= W; ++w)
        {
            AdeltaDot[w].clear();
            AdeltaGlobal[w].clear();
        }
    }
};

struct EdgeLocalCross : EdgeLocal
{
    std::unique_ptr<SCORETYPE *[]> Evals, Fvals, Gvals;
    std::unique_ptr<SCORETYPE ***[]> Esourcess, Fsourcess, Gsourcess;
    std::unique_ptr<SOURCESIZE *[]> Es_szs, Fs_szs, Gs_szs;
    std::unique_ptr<IDTYPE *[]> Eids, Fids, Gids;

    EdgeLocalCross(EdgeLocal &edge)
        : EdgeLocal(edge)
    {
    }

    SIZETYPE get_trn(SHORTSIZE ref_sz)
    {
        return 3 * (ref_sz + 1);
    }

    SIZETYPE get_tnn(QUERYSIZE Omax, SHORTSIZE ref_sz)
    {
        return (Omax + 1) * get_trn(ref_sz);
    }

    SIZETYPE get_tsn(QUERYSIZE Omax, SHORTSIZE ref_sz)
    {
        return SIZETYPE(Omax + 1) * 8 * (ref_sz + 1);
    }

    int64_t *get_Ss(SHORTSIZE ref_sz)
    {
        return new int64_t[3]{ref_sz, ref_sz, ref_sz};
    }

    SIZETYPE *get_SE(QUERYSIZE Omax, SHORTSIZE ref_sz)
    {
        SIZETYPE mid = SIZETYPE(Omax + 1) * 2 * (ref_sz + 1);
        return new SIZETYPE[3]{0, mid, get_tnn(Omax, ref_sz)};
    }

    SIZETYPE *get_steps()
    {
        return new SIZETYPE[2]{2, 4};
    }

    void apply_memory(QUERYSIZE Omax, SCORETYPE *&fpval, SCORETYPE **&fpsource, SCORETYPE ***&fpsources, SOURCESIZE *&fps_sz, IDTYPE *&fpid, SHORTSIZE *&fps, DOTTYPE *&fpn, SHORTSIZE ref_sz)
    {
        std::unique_ptr<int64_t []> Ss(get_Ss(ref_sz));
        extend_ss_ns(3, Ss.get(), n, fps, fpn);
        type_initial(fpval, {}, {&Evals, &Fvals, &Gvals}, Ss.get(), Omax);
        SCORETYPE ***fpsources_old = fpsources;
        type_initial(fpsources, {}, {&Esourcess, &Fsourcess, &Gsourcess}, Ss.get(), Omax);
        std::unique_ptr<SIZETYPE []> SE(get_SE(Omax, ref_sz));
        std::unique_ptr<SIZETYPE []> steps(get_steps());
        for (SIZETYPE i = 0; i < 2; ++i)
            for (SIZETYPE j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {}, {&Es_szs, &Fs_szs, &Gs_szs}, Ss.get(), Omax);
        type_initial(fpid, {}, {&Eids, &Fids, &Gids}, Ss.get(), Omax);
    }
};

struct EdgeLocalCircuit : EdgeLocal
{
    std::unique_ptr<SCORETYPE *[]> Evals, F0vals, G0vals, Gvals, D0vals, DXvals;
    SCORETYPE *Dvals;
    std::unique_ptr<SCORETYPE ***[]> Esourcess, F0sourcess, G0sourcess, Gsourcess;
    SCORETYPE ***Dsourcess;
    std::unique_ptr<SOURCESIZE *[]> Es_szs, F0s_szs, G0s_szs, Gs_szs;
    SOURCESIZE *Ds_szs;
    std::unique_ptr<IDTYPE *[]> Eids, F0ids, G0ids, Gids, D0ids, DXids;
    IDTYPE *Dids;
    std::unique_ptr<uint8_t []> DXbits;

    EdgeLocalCircuit(EdgeLocal &edge)
        : EdgeLocal(edge)
    {
    }

    SIZETYPE get_trn(SHORTSIZE ref_sz)
    {
        return 1 + 4 * (ref_sz + 1) + 2 * (tail->scc_sz - 1);
    }

    SIZETYPE get_tnn(QUERYSIZE Omax, SHORTSIZE ref_sz)
    {
        return (Omax + 1) * get_trn(ref_sz);
    }

    SIZETYPE get_tsn(QUERYSIZE Omax, SHORTSIZE ref_sz)
    {
        return SIZETYPE(Omax + 1) * (1 + 10 * (ref_sz + 1) + 3 * (tail->scc_sz - 1));
    }

    int64_t *get_Ss(SHORTSIZE ref_sz)
    {
        return new int64_t[6]{ref_sz, ref_sz, ref_sz, ref_sz, tail->scc_sz - 2, tail->scc_sz - 2};
    }

    SIZETYPE *get_SE(QUERYSIZE Omax, SHORTSIZE ref_sz)
    {
        return new SIZETYPE[6]{0, Omax + 1, SIZETYPE(Omax + 1) * (1 + 2 * (ref_sz + 1)), SIZETYPE(Omax + 1) * (1 + 4 * (ref_sz + 1))};
    }

    SIZETYPE *get_steps()
    {
        return new SIZETYPE[5]{1, 2, 3, 1, 2};
    }

    void apply_memory(QUERYSIZE Omax, SCORETYPE *&fpval, SCORETYPE **&fpsource, SCORETYPE ***&fpsources, SOURCESIZE *&fps_sz, IDTYPE *&fpid, SHORTSIZE *&fps, DOTTYPE *&fpn, SHORTSIZE ref_sz)
    {
        DXbits.reset(new uint8_t[(tail->scc_sz-1)*(Omax+1)]);
        std::unique_ptr<int64_t []> Ss(get_Ss(ref_sz));
        extend_ss_ns(1, n, fps, fpn);
        extend_ss_ns(6, Ss.get(), n, fps, fpn);
        type_initial(fpval, {&Dvals}, {&Evals, &F0vals, &G0vals, &Gvals, &D0vals, &DXvals}, Ss.get(), Omax);
        SCORETYPE ***fpsources_old = fpsources;
        type_initial(fpsources, {&Dsourcess}, {&Esourcess, &F0sourcess, &G0sourcess, &Gsourcess}, Ss.get(), Omax);
        fpsources += (tail->scc_sz - 1) * 2 * (Omax + 1);
        std::unique_ptr<SIZETYPE []> SE(get_SE(Omax, ref_sz));
        std::unique_ptr<SIZETYPE []> steps(get_steps());
        for (SIZETYPE i = 0; i < 3; ++i)
            for (SIZETYPE j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {&Ds_szs}, {&Es_szs, &F0s_szs, &G0s_szs, &Gs_szs}, Ss.get(), Omax);
        fps_sz += (tail->scc_sz - 1) * 2 * (Omax + 1);
        type_initial(fpid, {&Dids}, {&Eids, &F0ids, &G0ids, &Gids, &D0ids, &DXids}, Ss.get(), Omax);
    }
};

struct SNC
{
    SCORETYPE E;
    SCORETYPE F0;
    SCORETYPE G0;
    SCORETYPE hatG;
    NUCTYPE c;
    NUCTYPE idl;
    QUERYSIZE lambda;
    SIZETYPE sr1;
    SIZETYPE sr2;
    IDTYPE cid; // 0 means not try to add child yet
    SNC *itp;
    SNC *jump;

    SNC(SCORETYPE E_, SCORETYPE F0_, SCORETYPE G0_, SCORETYPE hatG_, NUCTYPE c_, NUCTYPE idl_, QUERYSIZE lambda_, SIZETYPE sr1_, SIZETYPE sr2_, SIZETYPE cid_, SNC *itp_, SNC *jump_)
        : E(E_), F0(F0_), G0(G0_), hatG(hatG_), c(c_), idl(idl_), lambda(lambda_), sr1(sr1_), sr2(sr2_), cid(cid_), itp(itp_), jump(jump_)
    {
    }
};

struct EdgeGlobalCircuit : Edge
{
    std::deque<SNC> sncs;
    std::unique_ptr<SCORETYPE *[]> D0vals;
    std::unique_ptr<IDTYPE *[]> D0ids;

    EdgeGlobalCircuit(Edge &edge)
        : Edge(edge)
    {
    }

    SIZETYPE get_trn()
    {
        return tail->scc_sz - 1;
    }

    SIZETYPE get_tnn(QUERYSIZE Omax)
    {
        return (Omax + 1) * get_trn();
    }

    int64_t *get_Ss()
    {
        return new int64_t[1]{tail->scc_sz - 2};
    }

    void apply_memory(SHORTSIZE Omax, SCORETYPE *&fpval, IDTYPE *&fpid, SHORTSIZE *&fps, DOTTYPE *&fpn)
    {
        std::unique_ptr<int64_t []> Ss(get_Ss());
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
    std::deque<Node> nodes;
    std::deque<EdgeLocalCross> local_crosses; // cannot be changed to vector, because emplace_back of vector may change the addresses of elements
    std::deque<EdgeLocalCircuit> local_circuits;
    std::deque<Edge> global_crosses;
    std::deque<EdgeGlobalCircuit> global_circuits;
    std::deque<Edge *> edges;

    std::deque<SCC> sccs;

    std::unique_ptr<SCORETYPE []> vals;
    std::unique_ptr<SCORETYPE *[]> sources;
    std::unique_ptr<SCORETYPE **[]> sourcess;
    std::unique_ptr<SOURCESIZE []> s_szs;
    std::unique_ptr<IDTYPE []> ids;
    std::unique_ptr<SHORTSIZE []> ss;
    std::unique_ptr<DOTTYPE []> ns;
    SIZETYPE trn, tnn, tsn;

    SCORETYPE *pQval;
    SCORETYPE ***pQsources;
    SOURCESIZE *pQs_sz;
    IDTYPE *pQid;

    boost::program_options::variables_map &vm;
    std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short;
    std::map<std::string, RankVec> &file2rankvec;
    std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SIZETYPE>> &file2long;
    std::map<std::string, std::ifstream> file2SA;

    void get_affine(std::vector<SCORETYPE> &g, SHORTSIZE seqlen, SCORETYPE initial, SCORETYPE u, SCORETYPE v)
    {
        g.resize(seqlen + 1);
        g[0] = initial;
        for (SIZETYPE s = 0; s < seqlen; ++s)
            if (s == 0)
                g[s + 1] = g[s] + v;
            else
                g[s + 1] = g[s] + u;
    }

    Graph(boost::program_options::variables_map &vm_, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short_, std::map<std::string, RankVec> &file2rankvec_, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SIZETYPE>> &file2long_) : vm(vm_), file2short(file2short_), file2rankvec(file2rankvec_), file2long(file2long_)
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
                for (NUCTYPE a = 2; a < RankVec::sigma; ++a)
                    for (NUCTYPE b = 2; b < RankVec::sigma; ++b)
                        if (a==b && a>2)
                            try{local.gamma[a][b]=std::stod(match);}catch(...){local.gamma[a][b]=1;}
                        else
                            try{local.gamma[a][b]=std::stod(mismatch);}catch(...){local.gamma[a][b]=-3;}
                local.name = name;
                for (Node &node : nodes)
                {
                    if (node.name == head)
                        local.head = &node;
                    if (node.name == tail)
                        local.tail = &node;
                }
                try{local.ve = std::stod(v);}catch(...){local.ve = -5;}
                local.vf = local.ve;
                try{local.ue = std::stod(u);}catch(...){local.ue = -2;}
                local.uf = local.ue;
                try{local.T = std::stod(T);}catch(...){local.T = -10;}
                try{local.min_score = std::stod(min_score);}catch(...){local.min_score = 20;}
                try{local.vfp = std::stod(vb);}catch(...){local.vfp = 0;}
                local.vfm = local.vfp;
                try{local.ufp = std::stod(ub);}catch(...){local.ufp = 0;}
                local.ufm = local.ufp;
                SHORTSIZE ref_sz = file2short[local.name].second;
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
                for (NUCTYPE a = 2; a < RankVec::sigma; ++a)
                    for (NUCTYPE b = 2; b < RankVec::sigma; ++b)
                        if (a==b && a>2)
                            try{global.gamma[a][b]=std::stod(match);}catch(...){global.gamma[a][b]=1;}
                        else
                            try{global.gamma[a][b]=std::stod(mismatch);}catch(...){global.gamma[a][b]=-3;}
                global.name = name;
                for (Node &node : nodes)
                {
                    if (node.name == head)
                        global.head = &node;
                    if (node.name == tail)
                        global.tail = &node;
                }
                try{global.ve = std::stod(v);}catch(...){global.ve = -5;}
                global.vf = global.ve;
                try{global.ue = std::stod(u);}catch(...){global.ue = -2;}
                global.uf = global.ue;
                try{global.T = std::stod(T);}catch(...){global.T = -10;}
                try{global.min_score = std::stod(min_score);}catch(...){global.min_score = 20;}

                globals.push_back(global);
            }

        std::vector<bool> visit(locals.size() + globals.size());
        RemoveEdge(nodes, false, visit, locals, globals);
        RemoveEdge(nodes, true, visit, locals, globals);
        SCCTS(locals, globals);
        for (SIZETYPE i = 0; i < sccs.size(); ++i)
            for (auto node : sccs[i].nodes)
            {
                node->scc_id = i;
                node->scc_sz = sccs[i].nodes.size();
            }
        edges.resize(locals.size() + globals.size());
        for (EdgeLocal &edge : locals)
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
        for (Edge &edge : globals)
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

    void DeepFirstLine(Node *node, bool reverse, std::vector<bool> &visit, std::list<EdgeLocal> &locals, std::list<Edge> &globals)
    {
        for (EdgeLocal &edge : locals)
            if (reverse ? edge.head == node : edge.tail == node)
                if (!visit[edge.n])
                {
                    visit[edge.n] = true;
                    DeepFirstLine(reverse ? edge.tail : edge.head, reverse, visit, locals, globals);
                }
        for (Edge &edge : globals)
            if (reverse ? edge.head == node : edge.tail == node)
                if (!visit[edge.n])
                {
                    visit[edge.n] = true;
                    DeepFirstLine(reverse ? edge.tail : edge.head, reverse, visit, locals, globals);
                }
    }

    void RemoveEdge(std::deque<Node> &nodes, bool reverse, std::vector<bool> &visit, std::list<EdgeLocal> &locals, std::list<Edge> &globals)
    {
        for (EdgeLocal &edge : locals)
            visit[edge.n] = false;
        for (Edge &edge : globals)
            visit[edge.n] = false;
        for (Node &node : nodes)
            if ((!reverse && node.is_root) || (reverse && node.is_target))
                DeepFirstLine(&node, reverse, visit, locals, globals);
        for (std::list<EdgeLocal>::iterator iter = locals.begin(); iter != locals.end();)
            if (!visit[iter->n])
                iter = locals.erase(iter);
            else
                ++iter;
        for (std::list<Edge>::iterator iter = globals.begin(); iter != globals.end();)
            if (!visit[iter->n])
                iter = globals.erase(iter);
            else
                ++iter;
    }

    void DeepFirst(Node *node, std::stack<Node *> &stack, SIZETYPE &id, std::vector<bool> &visit, std::vector<bool> &in_stack, std::vector<SIZETYPE> &disc, std::vector<SIZETYPE> &low, std::list<EdgeLocal> &locals, std::list<Edge> &globals)
    {
        DOTTYPE n = node->n;
        visit[n] = true;
        disc[n] = id;
        low[n] = id;
        ++id;
        stack.push(node);
        in_stack[n] = true;
        for (EdgeLocal &edge : locals)
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
        for (Edge &edge : globals)
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
        std::vector<SIZETYPE> disc(nodes.size()), low(nodes.size());
        for (SIZETYPE n = 0; n < nodes.size(); ++n)
        {
            visit[n] = false;
            in_stack[n] = false;
        }
        std::stack<Node *> stack;
        SIZETYPE id = 0;
        for (SIZETYPE n = 0; n < nodes.size(); ++n)
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
