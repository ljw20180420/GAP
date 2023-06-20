#ifndef TRACK_H
#define TRACK_H

#include "Align.h"

struct Track : Graph
{
    struct Peak
    {
        int n;
        std::string seqname;
        std::deque<std::pair<int64_t, int64_t>> sranges;
        std::pair<int64_t, int64_t> wrange;
        std::pair<double, double> valrange;
    };

    struct MegaRange
    {
        int fs;
        int s_sz;
    };

    std::ofstream fout_align, fout_extract, fout_fail, fout_death;
    int max_extract, diff_thres, max_range, min_seg_num, max_seg_num;
    int64_t max_mega;
    std::vector<Dot> dots;
    MonoDeque<Dot *> sources;
    std::map<Dot *, MegaRange> megadots_head;
    std::map<Dot *, MegaRange> megadots_tail;
    std::deque<std::pair<Dot *, std::deque<Dot *>>> megasources;

    Track(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel, int max_extract_, int64_t max_mega_, int diff_thres_, int max_range_, int min_seg_num_, int max_seg_num_)
    : Graph(argc, argv, file2seq, file2browheel), fout_align("alg"), fout_extract("ext"), fout_fail("fail"), fout_death("death"), max_extract(max_extract_), max_mega(max_mega_), diff_thres(diff_thres_), max_range(max_range_), min_seg_num(min_seg_num_), max_seg_num(max_seg_num_)
    {
    }

    void ReadTrack(std::deque<std::string> mg_files, std::deque<std::pair<std::string, std::string>> &reads_total)
    {
        std::vector<std::ifstream> fins;
        for (auto &file : mg_files)
            fins.emplace_back(file, std::ios::binary);
        int64_t max_id = 0;
        std::vector<std::string> Onames(fins.size());
        std::vector<int64_t> fzs(fins.size());
        for (int64_t i = 0; i < fins.size(); ++i)
        {
            int64_t id;
            fins[i].seekg(0, fins[i].end);
            fzs[i] = fins[i].tellg();
            if (fzs[i] > sizeof(int64_t))
            {
                fins[i].seekg(-sizeof(int64_t), fins[i].end);
                fins[i].read((char *)&id, sizeof(int64_t));
                max_id = id > max_id ? id : max_id;
                fins[i].seekg(0, fins[i].beg);
                int Onsz;
                fins[i].read((char *)&Onsz, sizeof(Onsz));
                Onames[i].resize(Onsz);
                fins[i].read((char *)Onames[i].data(), sizeof(char) * Onsz);
            }
        }
        dots.resize(max_id + 1);
        int64_t seq_num = 0;
        for (std::pair<std::string, std::string> &pair : reads_total)
        {
            ++seq_num;
            for (int64_t f = 0; f < fins.size(); ++f)
            {
                if (Onames[f] == pair.first)
                {
                    int64_t id = -1;
                    int64_t uid = 0;
                    do
                    {
                        ++id;
                        fins[f].read((char *)&dots[id].n, sizeof(Dot::n));
                        fins[f].read((char *)&dots[id].s, sizeof(Dot::s));
                        fins[f].read((char *)&dots[id].w, sizeof(Dot::w));
                        fins[f].read((char *)&dots[id].val, sizeof(Dot::val));
                        fins[f].read((char *)&dots[id].s_sz, sizeof(Dot::s_sz));
                        fins[f].read((char *)&dots[id].lambda, sizeof(Dot::lambda));
                        dots[id].fs = sources.offset;
                        for (int i = 0; i < dots[id].s_sz; ++i)
                        {
                            int64_t idd;
                            fins[f].read((char *)&idd, sizeof(idd));
                            sources.emplace_back(&dots[idd]);
                            uid = idd > uid ? idd : uid;
                        }
                    } while (id < uid);

                    BackTrack(&dots[0], pair.first, pair.second);

                    sources.clear();
                    if (fzs[f] - fins[f].tellg() > sizeof(int64_t))
                    {
                        int Onsz;
                        fins[f].read((char *)&Onsz, sizeof(Onsz));
                        Onames[f].resize(Onsz);
                        fins[f].read((char *)Onames[f].data(), sizeof(char) * Onsz);
                    }
                    break;
                }
            }
        }
    }

    void BackTrack(Dot *pQ, std::string &Oname, std::string &O)
    {
        megasources.clear();
        megadots_head.clear();
        megadots_tail.clear();
        megasources.emplace_back(pQ, std::deque<Dot *>());

        for (int64_t i = 0; i < megasources.size(); ++i)
            source_until(i);
        MegaRange Qrange = megadots_tail[pQ];

        MegaRange *pmegarange;
        std::deque<std::deque<int64_t>> extracts;

        for (int64_t i = Qrange.fs, next_mega = 0; i < Qrange.fs + Qrange.s_sz && extracts.size() < max_extract; ++i)
        {
            std::deque<int64_t> extract;
            extract.push_front(i);
            while (next_mega++ < max_mega && next_megadot(extract))
            {
                if (extract.size() % 2 == 0)
                    pmegarange = &megadots_tail[megasources[extract[0]].first];
                else
                    pmegarange = &megadots_head[megasources[extract[0]].first];
                if (pmegarange->s_sz == 0)
                {
                    bool addflag = true;
                    if (min_seg_num > 0 && max_seg_num > 0 && (extract.size() < 2 * min_seg_num || extract.size() > 2 * max_seg_num))
                        addflag = false;
                    else
                    {
                        for (auto &extract_pre : extracts)
                            if (extract_dis(extract_pre, extract) < diff_thres)
                            {
                                addflag = false;
                                break;
                            }
                    }
                    if (addflag)
                    {
                        extracts.push_back(extract);
                        if (extracts.size() == max_extract)
                            break;
                    }
                }
            }
        }

        extracts_to_aligns(extracts, Oname, O);
        output_extracts(extracts, Oname);

        if (extracts.empty())
        {
            if (Qrange.s_sz > 0)
            {
                fout_fail << Oname << '\n'
                      << O << '\n';
            }
            else
            {
                fout_death << Oname << '\n'
                      << O << '\n';
            }
        }
    }

    void source_until(int64_t i)
    {
        std::pair<Dot *, std::deque<Dot *>> &megasource = megasources[i];
        bool to_tail = i > 0 && megasource.second.empty();
        std::pair<std::map<Dot *, MegaRange>::iterator, bool> insertpair;
        if (to_tail)
            insertpair = megadots_head.emplace(megasource.first, MegaRange());
        else
            insertpair = megadots_tail.emplace(megasource.first, MegaRange());
        if (!insertpair.second)
            return;
        insertpair.first->second.fs = megasources.size();

        std::deque<Dot *> path;
        std::queue<Dot *> visited;
        megasource.first->visit = true;
        visited.push(megasource.first);
        path.push_back(megasource.first);
        do
        {
            if (to_tail)
            {
                for (int i = 0; i < path[0]->s_sz; ++i)
                    if (sources[path[0]->fs + i]->n < 0)
                    {
                        megasources.emplace_back(path[0], path);
                        break;
                    }
            }
            else if (path[0]->n >= 0 && path[0] != megasource.first)
                megasources.emplace_back(path[0], std::deque<Dot *>());

        } while (next_dot(path, visited, to_tail));
        insertpair.first->second.s_sz = megasources.size() - insertpair.first->second.fs;

        while (!visited.empty())
        {
            visited.front()->visit = false;
            visited.pop();
        }
    }

    bool next_dot(std::deque<Dot *> &path, std::queue<Dot *> &visited, bool to_tail)
    {
        int is = 0;
        while (!path.empty())
        {
            for (int i = is; i < path[0]->s_sz; ++i)
            {
                Dot *pdot = sources[path[0]->fs + i];
                if (pdot->visit || to_tail == (pdot->n < 0 || path[0]->n < 0))
                    continue;
                pdot->visit = true;
                visited.push(pdot);
                path.push_front(pdot);
                return true;
            }
            if (path.size() > 1)
            {
                for (is = 0; is < path[1]->s_sz; ++is)
                    if (sources[path[1]->fs + is] == path[0])
                    {
                        ++is;
                        break;
                    }
            }
            path.pop_front();
        }
        return false;
    }

    bool next_megadot(std::deque<int64_t> &extract)
    {
        MegaRange *pmegarange;
        int is = 0;
        while (!extract.empty())
        {
            if (extract.size() % 2 == 0)
                pmegarange = &megadots_tail[megasources[extract[0]].first];
            else
                pmegarange = &megadots_head[megasources[extract[0]].first];
            for (int i = is; i < pmegarange->s_sz; ++i)
            {
                int64_t idx = pmegarange->fs + i;
                Dot *pdot = megasources[idx].first;
                if (!legal_circle(pdot, extract))
                    continue;
                extract.push_front(idx);
                return true;
            }
            if (extract.size() > 1)
            {
                if (extract.size() % 2 == 0)
                    pmegarange = &megadots_head[megasources[extract[1]].first];
                else
                    pmegarange = &megadots_tail[megasources[extract[1]].first];
                is = extract[0] - pmegarange->fs + 1;
            }
            extract.pop_front();
        }
        return false;
    }

    bool legal_circle(Dot *pdot, std::deque<int64_t> &extract)
    {
        if (extract.size() % 2 == 0)
            return true;
        double summin_score = edges[pdot->n]->min_score;
        double sumdiffval = megasources[extract[0]].first->val - pdot->val;
        Node *new_node = edges[pdot->n]->tail;
        if (new_node == edges[pdot->n]->head && sumdiffval < summin_score)
            return false;
        for (int i = 1; i < extract.size(); i += 2)
        {
            summin_score += edges[megasources[extract[i]].first->n]->min_score;
            sumdiffval += megasources[extract[i + 1]].first->val - megasources[extract[i]].first->val;
            if (new_node == edges[megasources[extract[i]].first->n]->head && sumdiffval < summin_score)
                return false;
        }
        return true;
    }

    int64_t extract_dis(std::deque<int64_t> extract1, std::deque<int64_t> extract2)
    {
        int row = extract1.size() / 2 + 1, col = extract2.size() / 2 + 1;
        std::unique_ptr<int64_t[]> PD(new int64_t[row * col]), Ed(new int64_t[col - 1]), Fd(new int64_t[row - 1]);
        PD[0] = 0;
        for (int64_t i = 1; i < row; ++i)
        {
            if (!edges[megasources[extract1[2 * i - 1]].first->n]->global)
                Fd[i - 1] = megasources[extract1[2 * i - 1]].first->s - megasources[extract1[2 * i - 2]].first->s;
            else
                Fd[i - 1] = megasources[extract1[2 * i - 1]].first->lambda;
            PD[i] = PD[i - 1] + Fd[i - 1];
        }
        for (int64_t j = 1; j < col; ++j)
        {
            if (!edges[megasources[extract2[2 * j - 1]].first->n]->global)
                Ed[j - 1] = megasources[extract2[2 * j - 1]].first->s - megasources[extract2[2 * j - 2]].first->s;
            else
                Ed[j - 1] = megasources[extract2[2 * j - 1]].first->lambda;
            PD[j * row] = PD[(j - 1) * row] + Ed[j - 1];
        }
        for (int64_t i = 1; i < row; ++i)
            for (int64_t j = 1; j < col; ++j)
            {
                int64_t F = PD[j * row + i - 1] + Fd[i - 1];
                int64_t E = PD[(j - 1) * row + i] + Ed[j - 1];
                PD[j * row + i] = std::min(E, F);
                if (megasources[extract1[2 * i - 1]].first->n == megasources[extract2[2 * j - 1]].first->n)
                {
                    int64_t segdiff;
                    if (!edges[megasources[extract1[2 * i - 1]].first->n]->global)
                    {
                        int64_t inter = std::min(megasources[extract1[2 * i - 1]].first->s, megasources[extract2[2 * j - 1]].first->s) - std::max(megasources[extract1[2 * i - 2]].first->s, megasources[extract2[2 * j - 2]].first->s);
                        segdiff = Fd[i - 1] + Ed[j - 1] - 2 * (inter > 0 ? inter : 0);
                    }
                    else
                    {
                        std::map<std::string, std::deque<std::pair<int64_t, int64_t>>> seqranges1 = get_seqranges(megasources[extract1[2 * i - 1]].first), seqranges2 = get_seqranges(megasources[extract2[2 * j - 1]].first);
                        segdiff = global_dis(seqranges1, seqranges2);
                    }
                    PD[j * row + i] = std::min(PD[j * row + i], PD[(j - 1) * row + i - 1] + segdiff);
                }
            }
        return PD[row * col - 1];
    }

    int64_t global_dis(std::map<std::string, std::deque<std::pair<int64_t, int64_t>>> &seqranges1, std::map<std::string, std::deque<std::pair<int64_t, int64_t>>> &seqranges2)
    {
        int64_t max_inter = 0;
        for (auto it1 = seqranges1.begin(), it2 = seqranges2.begin(); it1 != seqranges1.end() && it2 != seqranges2.end();)
        {
            int cmp = it1->first.compare(it2->first);
            if (cmp < 0)
                ++it1;
            else if (cmp > 0)
                ++it2;
            else
            {
                for (auto &range1 : it1->second)
                    for (auto &range2 : it2->second)
                    {
                        int64_t inter = std::min(range1.second, range2.second) - std::max(range1.first, range2.first);
                        max_inter = std::max(max_inter, inter);
                    }
                ++it1;
                ++it2;
            }
        }
        std::pair<int64_t, int64_t> &range1 = (seqranges1.begin())->second.front(), &range2 = (seqranges2.begin())->second.front();
        return range1.second - range1.first + range2.second - range2.first - 2 * max_inter;
    }

    void extracts_to_aligns(std::deque<std::deque<int64_t>> &extracts, std::string &Oname, std::string &O)
    {
        for (int i = 0; i < extracts.size(); ++i)
        {
            std::string align_query, align_mid, align_ref;
            for (int j = 0; j < extracts[i].size(); j += 2)
            {
                std::deque<Dot *> &path = megasources[extracts[i][j]].second;
                if (!align_query.empty())
                {
                    align_query.push_back(' ');
                    align_mid.push_back(' ');
                    align_ref.push_back(' ');
                }
                EdgeLocal *local = static_cast<EdgeLocal *>(edges[path[0]->n]);
                EdgeGlobal *global = static_cast<EdgeGlobal *>(edges[path[0]->n]);
                if (!edges[path[0]->n]->global)
                {
                    align_query.append(path[0]->s, '-');
                    align_mid.append(path[0]->s, ' ');
                    align_ref.append(local->pnameseq->seq, 0, path[0]->s);
                }

                for (int i = 1; i < path.size(); ++i)
                {
                    bool E = path[i]->w > path[i - 1]->w, localF = path[i]->s > path[i - 1]->s && !edges[path[0]->n]->global, globalF = edges[path[i]->n]->global && path[i]->lambda > path[i - 1]->lambda;

                    if(!E && !localF && !globalF)
                        continue;

                    if (E)
                        align_query.push_back(O[path[i - 1]->w]);
                    else
                        align_query.push_back('-');
                    
                    if (localF)
                        align_ref.push_back(static_cast<EdgeLocal *>(edges[path[i]->n])->pnameseq->seq[path[i]->s - 1]);
                    else if (globalF)
                        align_ref.push_back(global->pbrowheel->sequence[global->pbrowheel->sequence.size() - 1 - path[i]->s]);
                    else
                        align_ref.push_back('-');

                    if (std::tolower(align_ref.back()) == std::tolower(align_query.back()))
                        align_mid.push_back('|');
                    else
                        align_mid.push_back(' ');
                }

                if (!edges[path[0]->n]->global)
                {
                    int sz = local->pnameseq->seq.size() - path.back()->s;
                    align_query.append(sz, '-');
                    align_mid.append(sz, ' ');
                    align_ref.append(local->pnameseq->seq, path.back()->s, sz);
                }
            }

            fout_align << Oname << '\t' << i << '\n'
                       << align_query << '\n'
                       << align_mid << '\n'
                       << align_ref << '\n';
        }
    }

    void output_extracts(std::deque<std::deque<int64_t>> &extracts, std::string &Oname)
    {
        for (int i = 0; i < extracts.size(); ++i)
        {
            fout_extract << Oname << '\t' << i << '\n';
            for (int j = 0; j < extracts[i].size(); j += 2)
            {
                Dot *pdot1 = megasources[extracts[i][j]].first;
                Dot *pdot2 = megasources[extracts[i][j + 1]].first;
                Edge *edge = edges[pdot2->n];
                if (edge->global)
                {
                    std::map<std::string, std::deque<std::pair<int64_t, int64_t>>> seqranges = get_seqranges(pdot2);
                    fout_extract << edge->name << ':';
                    for (auto &seqrange : seqranges)
                        for (auto &range : seqrange.second)
                            fout_extract << seqrange.first << '\t' << range.first << '\t';
                    fout_extract << pdot1->w << '\t' << pdot1->val << '\n';
                    fout_extract << edge->name << ':';
                    for (auto &seqrange : seqranges)
                        for (auto &range : seqrange.second)
                            fout_extract << seqrange.first << '\t' << range.second << '\t';
                    fout_extract << pdot2->w << '\t' << pdot2->val << '\n';
                }
                else
                {
                    fout_extract << edge->name << '\t' << pdot1->s << '\t' << pdot1->w << '\t' << pdot1->val << '\n'
                                 << edge->name << '\t' << pdot2->s << '\t' << pdot2->w << '\t' << pdot2->val << '\n';
                }
            }
        }
    }

    std::map<std::string, std::deque<std::pair<int64_t, int64_t>>> get_seqranges(Dot *pdot)
    {
        std::map<std::string, std::deque<std::pair<int64_t, int64_t>>> seqranges;
        BroWheel *pbrowheel = static_cast<EdgeGlobal *>(edges[pdot->n])->pbrowheel;
        int64_t sr1 = 0, sr2 = pbrowheel->sequence.size() - 1;
        for (int64_t s = pbrowheel->sequence.size() - pdot->s + pdot->lambda - 2; s >= pbrowheel->sequence.size() - 1 - pdot->s; --s)
            pbrowheel->PreRange(sr1, sr2, BroWheel::base2int(pbrowheel->sequence[s]));
        for (int64_t sx = sr1; sx <= sr2 && sx - sr1 < max_range; ++sx)
        {
            int64_t s = pbrowheel->SimSuffix(sx);
            int64_t l = 0, r = pbrowheel->name_cumlen.size();
            while (r - l > 1)
            {
                int64_t m = (l + r) / 2;
                if (pbrowheel->name_cumlen[m].second > s)
                    r = m;
                else
                    l = m;
            }
            if (r < pbrowheel->name_cumlen.size())
                s = pbrowheel->name_cumlen[r].second - 1 - s;
            else
                s = pbrowheel->sequence.size() - 1 - s;
            seqranges[pbrowheel->name_cumlen[l].first].emplace_back(s - pdot->lambda, s);
        }
        return seqranges;
    }
};

#endif