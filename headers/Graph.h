#ifndef GRAPH_H
#define GRAPH_H

#include <stack>
#include <cfloat>
#include <list>
#include "BroWheel.h"

constexpr static const double inf = std::numeric_limits<double>::infinity();

struct Dot
{
    static const int DotQ = -3, DotAbar = -2, DotB = -1;

    bool visit;
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
    bool global;
    std::string name;
    double gamma[sigma + 1][sigma + 1] = {
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -inf, -inf, -inf, -inf, -inf},
        {-inf, -inf, -3, -3, -3, -3, -3},
        {-inf, -inf, -3, 1, -3, -3, -3},
        {-inf, -inf, -3, -3, 1, -3, -3},
        {-inf, -inf, -3, -3, -3, 1, -3},
        {-inf, -inf, -3, -3, -3, -3, 1}};
    double ve, ue, vf, uf, T, dT, min_score;
    Node *tail;
    Node *head;
    std::vector<double> gf;
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
    NameSeq *pnameseq;
    double vfp;
    double ufp;
    double vfm;
    double ufm;
    double gfm;
    std::vector<double> gfpT;

    EdgeLocal(Edge &edge)
        : Edge(edge)
    {
    }
};

struct TrackTree
{
    enum enumEFG
    {
        enumE,
        enumF,
        enumG
    };

    int W;
    std::deque<Dot> dots;
    int64_t idx, pidx;
    int shiftE, shiftF, shiftG, shiftPE, shiftPF, shiftPG;

    void setidx(int64_t idx_)
    {
        if (idx_ == 0)
            pidx = -1;
        else
        {
            pidx = idx;
            shiftPE = shiftE;
            shiftPF = shiftF;
            shiftPG = shiftG;
        }
        idx = idx_;
        shiftE = 3 * idx * (W + 1);
        shiftF = shiftE + (W + 1);
        shiftG = shiftF + (W + 1);
    }
    struct TrackNode
    {
        int64_t cidxs[5];
        int tau = 0;

        TrackNode(int64_t cidx_)
        {
            for (int i = 0; i < 5; ++i)
                cidxs[i] = cidx_;
        }
    };

    std::deque<TrackNode> tracknodes;

    void emplace_back(int64_t cidx_, int n_, int64_t s_, int lambda_)
    {
        tracknodes.emplace_back(cidx_);
        dots.resize(tracknodes.size() * (enumG + 1) * (W + 1));
        for (int k = (tracknodes.size() - 1) * (enumG + 1) * (W + 1), t = enumE; t <= enumG; ++t)
            for (int w = 0; w <= W; ++w, ++k)
            {
                dots[k].n = n_;
                dots[k].s = s_;
                dots[k].w = w;
                dots[k].lambda = lambda_;
            }
    }

    void clear()
    {
        dots.clear();
        tracknodes.clear();
    }
};

struct EdgeGlobal : Edge
{
    bool reverse_complement;
    BroWheel *pbrowheel;
    TrackTree tracktree;

    EdgeGlobal(Edge &edge)
        : Edge(edge)
    {
    }
};

template <typename T>
void type_initial(T *&Ts, std::initializer_list<T **> pTss, std::initializer_list<std::unique_ptr<T *[]> *> pTsss, int *Ss, int Omax)
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
struct EdgeGlobalCross;
struct EdgeGlobalCircuit;

struct Node
{
    bool is_root, is_target;
    int n, scc_id, scc_sz;
    std::string name;
    double ve, ue;
    std::vector<EdgeLocalCross *> out_local_crosses, in_local_crosses;
    std::vector<EdgeLocalCircuit *> out_local_circuits, in_local_circuits;
    std::vector<EdgeGlobalCross *> out_global_crosses, in_global_crosses;
    std::vector<EdgeGlobalCircuit *> out_global_circuits, in_global_circuits;

    std::unique_ptr<bool *[]> Avisits;
    bool *Bvisits, *pAbarvisit;
    std::unique_ptr<double *[]> Avals;
    double *Bvals, *pAbarval;
    std::unique_ptr<double ***[]> Asourcess;
    double ***Bsourcess, ***pAbarsources;
    std::unique_ptr<int *[]> As_szs;
    int *Bs_szs, *pAbars_sz;
    std::unique_ptr<int64_t *[]> Aids;
    int64_t *Bids, *pAbarid;

    std::unique_ptr<std::deque<double *>[]> AdeltaDot;

    struct GlobalSuffix
    {
        EdgeGlobal *edge;
        int64_t start;
        int lambda;

        GlobalSuffix(EdgeGlobal *edge_, int64_t start_, int lambda_)
            : edge(edge_), start(start_), lambda(lambda_)
        {
        }
    };

    std::unique_ptr<std::deque<GlobalSuffix>[]> AdeltaGlobal;

    Node(std::string name_, double ve_, double ue_, int n_)
        : name(name_), ve(ve_), ue(ue_), n(n_), is_root(false), is_target(false)
    {
    }

    int get_trn()
    {
        return 1 + scc_sz;
    }

    int64_t get_tnn(int Omax)
    {
        return int64_t(Omax + 1) * get_trn();
    }

    int64_t get_tsn(int Omax)
    {
        int Aso = 1 + 2 * in_local_circuits.size() + in_global_circuits.size();
        return int64_t(Omax + 1) * (2 + Aso * (scc_sz - 1));
    }

    int *get_Ss()
    {
        return new int[1]{scc_sz - 1};
    }

    int64_t *get_SE(int Omax)
    {
        return new int64_t[4]{0, Omax + 1, 2 * (Omax + 1), get_tnn(Omax)};
    }

    int *get_steps()
    {
        int Aso = 1 + 2 * in_local_circuits.size() + in_global_circuits.size();
        return new int[3]{2, 0, Aso};
    }

    void apply_memory(int Omax, bool *&fpvisit, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn, bool *&rpvisit, double *&rpval, double ***&rpsources, int *&rps_sz, int64_t *&rpid)
    {
        std::unique_ptr<int[]> Ss(new int[1]{scc_sz - 1});
        extend_ss_ns(1, -1, fps, fpn);
        extend_ss_ns(1, Ss.get(), Dot::nidx_trans(n), fps, fpn);
        type_initial(fpvisit, {&Bvisits}, {&Avisits}, Ss.get(), Omax);
        pAbarvisit = --rpvisit;
        type_initial(fpval, {&Bvals}, {&Avals}, Ss.get(), Omax);
        pAbarval = --rpval;
        if (is_root)
            *pAbarval = 0;
        else
            *pAbarval = -inf;
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {&Bsourcess}, {&Asourcess}, Ss.get(), Omax);
        pAbarsources = --rpsources;
        std::unique_ptr<int64_t[]> SE(get_SE(Omax));
        std::unique_ptr<int[]> steps(get_steps());
        for (int i = 0; i < 3; ++i)
            for (int64_t j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {&Bs_szs}, {&As_szs}, Ss.get(), Omax);
        pAbars_sz = --rps_sz;
        *pAbars_sz = 0;
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

    void updateA0(int w, double source_val, EdgeGlobal *edge_, int64_t start_, int lambda_)
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
    std::unique_ptr<bool *[]> Evisits, Fvisits, Gvisits;
    std::unique_ptr<double *[]> Evals, Fvals, Gvals;
    std::unique_ptr<double ***[]> Esourcess, Fsourcess, Gsourcess;
    std::unique_ptr<int *[]> Es_szs, Fs_szs, Gs_szs;
    std::unique_ptr<int64_t *[]> Eids, Fids, Gids;

    EdgeLocalCross(EdgeLocal &edge)
        : EdgeLocal(edge)
    {
    }

    int get_trn()
    {
        return 3 * (pnameseq->seq.size() + 1);
    }

    int64_t get_tnn(int Omax)
    {
        return int64_t(Omax + 1) * get_trn();
    }

    int64_t get_tsn(int Omax)
    {
        return int64_t(Omax + 1) * 8 * (pnameseq->seq.size() + 1);
    }

    int *get_Ss()
    {
        int S = pnameseq->seq.size();
        return new int[3]{S, S, S};
    }

    int64_t *get_SE(int Omax)
    {
        int64_t mid = int64_t(Omax + 1) * 2 * (pnameseq->seq.size() + 1);
        return new int64_t[3]{0, mid, get_tnn(Omax)};
    }

    int *get_steps()
    {
        return new int[2]{2, 4};
    }

    void apply_memory(int Omax, bool *&fpvisit, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn)
    {
        std::unique_ptr<int[]> Ss(get_Ss());
        extend_ss_ns(3, Ss.get(), n, fps, fpn);
        type_initial(fpvisit, {}, {&Evisits, &Fvisits, &Gvisits}, Ss.get(), Omax);
        type_initial(fpval, {}, {&Evals, &Fvals, &Gvals}, Ss.get(), Omax);
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {}, {&Esourcess, &Fsourcess, &Gsourcess}, Ss.get(), Omax);
        std::unique_ptr<int64_t[]> SE(get_SE(Omax));
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
    std::unique_ptr<bool *[]> Evisits, F0visits, G0visits, Gvisits, D0visits, DXvisits;
    bool *Dvisits;
    std::unique_ptr<double *[]> Evals, F0vals, G0vals, Gvals, D0vals, DXvals;
    double *Dvals;
    std::unique_ptr<double ***[]> Esourcess, F0sourcess, G0sourcess, Gsourcess, D0sourcess, DXsourcess;
    double ***Dsourcess;
    std::unique_ptr<int *[]> Es_szs, F0s_szs, G0s_szs, Gs_szs, D0s_szs, DXs_szs;
    int *Ds_szs;
    std::unique_ptr<int64_t *[]> Eids, F0ids, G0ids, Gids, D0ids, DXids;
    int64_t *Dids;

    EdgeLocalCircuit(EdgeLocal &edge)
        : EdgeLocal(edge)
    {
    }

    int get_trn()
    {
        return 1 + 4 * (pnameseq->seq.size() + 1) + 2 * (tail->scc_sz - 1);
    }

    int64_t get_tnn(int Omax)
    {
        return int64_t(Omax + 1) * get_trn();
    }

    int64_t get_tsn(int Omax)
    {
        return int64_t(Omax + 1) * (1 + 10 * (pnameseq->seq.size() + 1) + 3 * (tail->scc_sz - 1));
    }

    int *get_Ss()
    {
        int S = pnameseq->seq.size();
        int scc_sz = tail->scc_sz;
        return new int[6]{S, S, S, S, scc_sz - 2, scc_sz - 2};
    }

    int64_t *get_SE(int Omax)
    {
        int S = pnameseq->seq.size();
        int scc_sz = tail->scc_sz;
        return new int64_t[6]{0, Omax + 1, int64_t(Omax + 1) * (1 + 2 * (S + 1)), int64_t(Omax + 1) * (1 + 4 * (S + 1)), int64_t(Omax + 1) * (1 + 4 * (S + 1) + (scc_sz - 1)), int64_t(Omax + 1) * (1 + 4 * (S + 1) + 2 * (scc_sz - 1))};
    }

    int *get_steps()
    {
        return new int[5]{1, 2, 3, 1, 2};
    }

    void apply_memory(int Omax, bool *&fpvisit, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn)
    {
        std::unique_ptr<int[]> Ss(get_Ss());
        extend_ss_ns(1, n, fps, fpn);
        extend_ss_ns(6, Ss.get(), n, fps, fpn);
        type_initial(fpvisit, {&Dvisits}, {&Evisits, &F0visits, &G0visits, &Gvisits, &D0visits, &DXvisits}, Ss.get(), Omax);
        type_initial(fpval, {&Dvals}, {&Evals, &F0vals, &G0vals, &Gvals, &D0vals, &DXvals}, Ss.get(), Omax);
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {&Dsourcess}, {&Esourcess, &F0sourcess, &G0sourcess, &Gsourcess, &D0sourcess, &DXsourcess}, Ss.get(), Omax);
        std::unique_ptr<int64_t[]> SE(get_SE(Omax));
        std::unique_ptr<int[]> steps(get_steps());
        for (int i = 0; i < 5; ++i)
            for (int64_t j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {&Ds_szs}, {&Es_szs, &F0s_szs, &G0s_szs, &Gs_szs, &D0s_szs, &DXs_szs}, Ss.get(), Omax);
        type_initial(fpid, {&Dids}, {&Eids, &F0ids, &G0ids, &Gids, &D0ids, &DXids}, Ss.get(), Omax);
    }
};

struct EdgeGlobalCross : EdgeGlobal
{
    EdgeGlobalCross(EdgeGlobal &edge)
        : EdgeGlobal(edge)
    {
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

struct EdgeGlobalCircuit : EdgeGlobal
{
    std::deque<SNC> sncs;
    std::unique_ptr<bool *[]> D0visits;
    std::unique_ptr<double *[]> D0vals;
    std::unique_ptr<double ***[]> D0sourcess;
    std::unique_ptr<int *[]> D0s_szs;
    std::unique_ptr<int64_t *[]> D0ids;

    EdgeGlobalCircuit(EdgeGlobal &edge)
        : EdgeGlobal(edge)
    {
    }

    int get_trn()
    {
        return tail->scc_sz - 1;
    }

    int64_t get_tnn(int Omax)
    {
        return int64_t(Omax + 1) * get_trn();
    }

    int64_t get_tsn(int Omax)
    {
        return int64_t(Omax + 1) * (tail->scc_sz - 1);
    }

    int *get_Ss()
    {
        return new int[1]{tail->scc_sz - 2};
    }

    int64_t *get_SE(int Omax)
    {
        return new int64_t[2]{0, get_tnn(Omax)};
    }

    int *get_steps()
    {
        return new int[1]{1};
    }

    void apply_memory(int Omax, bool *&fpvisit, double *&fpval, double **&fpsource, double ***&fpsources, int *&fps_sz, int64_t *&fpid, int *&fps, int *&fpn)
    {
        std::unique_ptr<int[]> Ss(get_Ss());
        extend_ss_ns(1, Ss.get(), n, fps, fpn);
        type_initial(fpvisit, {}, {&D0visits}, Ss.get(), Omax);
        type_initial(fpval, {}, {&D0vals}, Ss.get(), Omax);
        double ***fpsources_old = fpsources;
        type_initial(fpsources, {}, {&D0sourcess}, Ss.get(), Omax);
        std::unique_ptr<int64_t[]> SE(get_SE(Omax));
        std::unique_ptr<int[]> steps(get_steps());
        for (int i = 0; i < 1; ++i)
            for (int64_t j = SE[i]; j < SE[i + 1]; ++j, fpsource += steps[i])
                fpsources_old[j] = fpsource;
        type_initial(fps_sz, {}, {&D0s_szs}, Ss.get(), Omax);
        type_initial(fpid, {}, {&D0ids}, Ss.get(), Omax);
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

std::string trim_whitespace(std::string str)
{
    size_t i = str.find_first_not_of(" \n\r\t\f\v"), j = str.find_last_not_of(" \n\r\t\f\v");
    return str.substr(i, j - i + 1);
}

std::deque<std::string> split_string(std::string str, char delimiters, char open_range, char close_range)
{
    std::deque<std::string> substrs;
    int depth = 0, prei = -1;
    for (int i = 0; i < str.size(); ++i)
    {
        if (str[i] == delimiters)
        {
            if (depth == 0)
            {
                int len = i - prei - 1;
                if (len > 0)
                    substrs.emplace_back(trim_whitespace(str.substr(prei + 1, len)));
                prei = i;
            }
        }
        else if (str[i] == open_range)
            ++depth;
        else if (str[i] == close_range)
            --depth;
    }
    return substrs;
}

std::string remove_open_close(std::string str, char open_range, char close_range)
{
    size_t i = str.find_first_of(open_range), j = str.find_last_of(close_range);
    return str.substr(i + 1, j - i - 1);
}

std::map<std::string, std::string> get_handles(std::deque<std::string> strs)
{
    std::map<std::string, std::string> handles_to_strs;
    for (std::string &str : strs)
        for (int i = 0; i < str.size(); ++i)
            if (str[i] == '=')
            {
                handles_to_strs[trim_whitespace(str.substr(0, i))] = trim_whitespace(str.substr(i + 1, str.size() - i - 1));
                break;
            }
    return handles_to_strs;
}

struct Graph
{
    const static int sigma = 6;

    std::deque<Node> nodes;
    std::vector<Node *> roots;
    std::vector<Node *> targets;
    std::deque<EdgeLocalCross> local_crosses; // cannot be changed to vector, because emplace_back of vector may change the addresses of elements
    std::deque<EdgeLocalCircuit> local_circuits;
    std::deque<EdgeGlobalCross> global_crosses;
    std::deque<EdgeGlobalCircuit> global_circuits;
    std::deque<Edge *> edges;

    std::deque<SCC> sccs;

    std::unique_ptr<bool[]> visits;
    std::unique_ptr<double[]> vals;
    std::unique_ptr<double *[]> sources;
    std::unique_ptr<double **[]> sourcess;
    std::unique_ptr<int[]> s_szs;
    std::unique_ptr<int64_t[]> ids;
    std::unique_ptr<int[]> ss, ns;
    int trn;
    int64_t tnn, tsn;

    bool *pQvisit;
    double *pQval;
    double ***pQsources;
    int *pQs_sz;
    int64_t *pQid;

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

    Graph(std::string argfile, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel)
    {
        std::string data;
        std::ifstream fin(argfile);
        std::getline(fin, data, '\0');
        std::deque<std::string> strs = split_string(data, ';', '{', '}');
        std::map<std::string, std::string> hts = get_handles(strs);

        std::deque<std::string> settings = split_string(remove_open_close(hts["nodes"], '{', '}'), ';', '{', '}');
        for (std::string &setting : settings)
        {
            std::map<std::string, std::string> htnp = get_handles(split_string(remove_open_close(setting, '{', '}'), ';', '\0', '\0'));
            nodes.emplace_back(htnp["name"], str2double(htnp["ve"]), str2double(htnp["ue"]), nodes.size());
        }

        settings = split_string(remove_open_close(hts["roots"], '{', '}'), ';', '\0', '\0');
        for (auto &node : nodes)
            if (std::find(settings.begin(), settings.end(), node.name) != settings.end())
            {
                roots.push_back(&node);
                node.is_root = true;
                if (roots.size() == settings.size())
                    break;
            }
        settings = split_string(remove_open_close(hts["targets"], '{', '}'), ';', '\0', '\0');
        for (auto &node : nodes)
            if (std::find(settings.begin(), settings.end(), node.name) != settings.end())
            {
                targets.push_back(&node);
                node.is_target = true;
                if (targets.size() == settings.size())
                    break;
            }

        std::list<EdgeLocal> locals;
        settings = split_string(remove_open_close(hts["locals"], '{', '}'), ';', '{', '}');
        for (std::string &setting : settings)
        {
            Edge edge;
            edge.n = locals.size();
            edge.global = false;
            std::map<std::string, std::string> htp = get_handles(split_string(remove_open_close(setting, '{', '}'), ';', '{', '}'));
            if (htp["gamma"] != "default")
            {
                std::deque<std::string> pair_scores = split_string(remove_open_close(htp["gamma"], '{', '}'), ';', '\0', '\0');
                for (int a = 0; a < 7; ++a)
                    for (int b = 0; b < 7; ++b)
                        edge.gamma[a][b] = str2double(pair_scores[7 * a + b]);
            }
            edge.name = htp["name"];
            edge.ve = str2double(htp["ve"]);
            edge.ue = str2double(htp["ue"]);
            edge.vf = str2double(htp["vf"]);
            edge.uf = str2double(htp["uf"]);
            edge.T = str2double(htp["T"]);
            edge.dT = str2double(htp["dT"]);
            edge.min_score = str2double(htp["min_score"]);
            get_affine(edge.gf, file2seq[edge.name].seq.size(), 0, edge.uf, edge.vf);
            for (auto &node : nodes)
            {
                if (node.name == htp["tail"])
                    edge.tail = &node;
                if (node.name == htp["head"])
                    edge.head = &node;
            }
            locals.emplace_back(edge);
            EdgeLocal &local = locals.back();
            local.pnameseq = &file2seq[edge.name];
            local.vfp = str2double(htp["vfp"]);
            local.ufp = str2double(htp["ufp"]);
            local.vfm = str2double(htp["vfm"]);
            local.ufm = str2double(htp["ufm"]);
            local.gfm = local.vfm + (local.pnameseq->seq.size() - 1) * local.ufm;
            get_affine(local.gfpT, local.pnameseq->seq.size(), local.T, local.ufp, local.vfp);
        }

        std::list<EdgeGlobal> globals;
        settings = split_string(remove_open_close(hts["globals"], '{', '}'), ';', '{', '}');
        for (std::string &setting : settings)
        {
            Edge edge;
            edge.n = locals.size() + globals.size();
            edge.global = true;
            std::map<std::string, std::string> htp = get_handles(split_string(remove_open_close(setting, '{', '}'), ';', '{', '}'));
            if (htp["gamma"] != "default")
            {
                std::deque<std::string> pair_scores = split_string(remove_open_close(htp["gamma"], '{', '}'), ';', '\0', '\0');
                for (int a = 0; a < 7; ++a)
                    for (int b = 0; b < 7; ++b)
                        edge.gamma[a][b] = str2double(pair_scores[7 * a + b]);
            }
            edge.name = htp["name"];
            edge.ve = str2double(htp["ve"]);
            edge.ue = str2double(htp["ue"]);
            edge.vf = str2double(htp["vf"]);
            edge.uf = str2double(htp["uf"]);
            edge.T = str2double(htp["T"]);
            edge.dT = str2double(htp["dT"]);
            edge.min_score = str2double(htp["min_score"]);
            get_affine(edge.gf, file2seq[edge.name].seq.size(), 0, edge.uf, edge.vf);
            for (auto &node : nodes)
            {
                if (node.name == htp["tail"])
                    edge.tail = &node;
                if (node.name == htp["head"])
                    edge.head = &node;
            }
            globals.emplace_back(edge);
            globals.back().pbrowheel = &file2browheel[edge.name];
            transform(htp["reverse_complement"].begin(), htp["reverse_complement"].end(), htp["reverse_complement"].begin(), ::tolower);
            if (htp["reverse_complement"] == "true")
                globals.back().reverse_complement = true;
            else if (htp["reverse_complement"] == "false")
                globals.back().reverse_complement = false;
            else
            {
                std::cerr << "reverse_complement must be true or false\n";
                exit(-1);
            }
        }
        std::vector<bool> visit(locals.size() + globals.size());
        RemoveEdge(roots, false, visit, locals, globals);
        RemoveEdge(targets, true, visit, locals, globals);
        SCCTS(locals, globals);
        for (int i = 0; i < sccs.size(); ++i)
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
                edge.tail->out_local_crosses.emplace_back(&local_crosses.back());
                edge.head->in_local_crosses.emplace_back(&local_crosses.back());
                edges[edge.n] = &local_crosses.back();
            }
            else
            {
                local_circuits.emplace_back(edge);
                sccs[edge.tail->scc_id].local_circuits.emplace_back(&local_circuits.back());
                edge.tail->out_local_circuits.emplace_back(&local_circuits.back());
                edge.head->in_local_circuits.emplace_back(&local_circuits.back());
                edges[edge.n] = &local_circuits.back();
            }
        for (auto &edge : globals)
            if (edge.tail->scc_id != edge.head->scc_id)
            {
                global_crosses.emplace_back(edge);
                sccs[edge.tail->scc_id].global_crosses.emplace_back(&global_crosses.back());
                edge.tail->out_global_crosses.emplace_back(&global_crosses.back());
                edge.head->in_global_crosses.emplace_back(&global_crosses.back());
                edges[edge.n] = &global_crosses.back();
            }
            else
            {
                global_circuits.emplace_back(edge);
                sccs[edge.tail->scc_id].global_circuits.emplace_back(&global_circuits.back());
                edge.tail->out_global_circuits.emplace_back(&global_circuits.back());
                edge.head->in_global_circuits.emplace_back(&global_circuits.back());
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
                iter = locals.erase(iter);
            else
                ++iter;
        for (auto iter = globals.begin(); iter != globals.end();)
            if (!visit[iter->n])
                iter = globals.erase(iter);
            else
                ++iter;
    }

    void DeepFirst(Node *node, std::stack<Node *> &stack, int &id, std::vector<bool> &visit, std::vector<bool> &in_stack, std::vector<int> &disc, std::vector<int> &low, std::list<EdgeLocal> &locals, std::list<EdgeGlobal> &globals)
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

    void SCCTS(std::list<EdgeLocal> &locals, std::list<EdgeGlobal> &globals)
    {
        std::vector<bool> visit(nodes.size()), in_stack(nodes.size());
        std::vector<int> disc(nodes.size()), low(nodes.size());
        for (int64_t n = 0; n < nodes.size(); ++n)
        {
            visit[n] = false;
            in_stack[n] = false;
        }
        std::stack<Node *> stack;
        int id = 0;
        for (int64_t n = 0; n < nodes.size(); ++n)
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
