#ifndef TRACK_H
#define TRACK_H

#include "Align.h"

struct Track : Graph
{
    struct MegaRange
    {
        uint64_t fs;
        uint64_t s_sz;
    };

    std::ofstream fout_align, fout_extract, fout_fail, fout_death;
    std::vector<Dot> dots;
    MonoDeque<Dot *> sources;
    std::map<Dot *, MegaRange> megadots_head;
    std::map<Dot *, MegaRange> megadots_tail;
    std::deque<std::pair<Dot *, std::deque<Dot *>>> megasources;
    std::map<std::string, std::vector<std::pair<std::string, uint64_t>>> file2cumlen;

    Track(boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2short, std::map<std::string, RankVec> &file2rankvec, std::map<std::string, std::pair<std::unique_ptr<uint8_t[]>, uint64_t>> &file2long)
    : Graph(vm, file2short, file2rankvec, file2long), fout_align("alg"), fout_extract("ext"), fout_fail("fail"), fout_death("death")
    {
        for (std::string file : vm["longs"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2cumlen.count(file))
                continue;

            std::string refname;
            uint64_t cumlen;
            for (std::ifstream fin(file+".cumlen"); fin >> refname >> cumlen;)
                file2cumlen[file].emplace_back(refname, cumlen);
        }
    }

    void ReadTrack(std::deque<std::string> mg_files)
    {
        std::vector<std::ifstream> fins;
        std::vector<uintmax_t> fzs;
        std::vector<std::string> Onames;
        decltype(Dot::id) max_id = 0, id;
        for (auto &mg_file : mg_files)
        {
            uintmax_t fz = std::filesystem::file_size(mg_file);
            if (fz <= sizeof(Dot::id))
                continue;
            fzs.push_back(fz);
            fins.emplace_back(mg_file, std::ios::binary);
            fins.back().seekg(-sizeof(Dot::id), fins.back().end);
            fins.back().read((char *)&id, sizeof(Dot::id));
            max_id = id > max_id ? id : max_id;
            fins.back().seekg(0, fins.back().beg);
            uint64_t Onsz;
            fins.back().read((char *)&Onsz, sizeof(Onsz));
            Onames.emplace_back(Onsz, 'x');
            fins.back().read((char *)Onames.back().data(), Onsz);
        }

        dots.resize(max_id + 1);
        std::ifstream finput(vm["input"].as<std::string>());
        for (std::string name, read; std::getline(std::getline(finput, name), read); )
        {
            for (uint64_t f = 0; f < fins.size(); ++f)
            {
                if (Onames[f] == name)
                {
                    for (uint64_t id = 0, uid = 0; id <= uid; ++id)
                    {
                        fins[f].read((char *)&dots[id].n, sizeof(Dot::n));
                        fins[f].read((char *)&dots[id].s, sizeof(Dot::s));
                        fins[f].read((char *)&dots[id].w, sizeof(Dot::w));
                        fins[f].read((char *)&dots[id].val, sizeof(Dot::val));
                        fins[f].read((char *)&dots[id].s_sz, sizeof(Dot::s_sz));
                        fins[f].read((char *)&dots[id].lambda, sizeof(Dot::lambda));
                        dots[id].fs = sources.offset;
                        for (uint64_t i = 0; i < dots[id].s_sz; ++i)
                        {
                            uint64_t idd;
                            fins[f].read((char *)&idd, sizeof(idd));
                            sources.emplace_back(&dots[idd]);
                            uid = idd > uid ? idd : uid;
                        }
                    }

                    BackTrack(&dots[0], name, read);
                    sources.clear();

                    if (fzs[f]-fins[f].tellg()>sizeof(Dot::id))
                    {
                        uint64_t Onsz;
                        fins[f].read((char *)&Onsz, sizeof(Onsz));
                        Onames[f].resize(Onsz);
                        fins[f].read((char *)Onames[f].data(), Onsz);
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

        for (uint64_t i = 0; i < megasources.size(); ++i)
            source_until(i);
        MegaRange Qrange = megadots_tail[pQ];

        MegaRange *pmegarange;
        std::deque<std::deque<int64_t>> extracts;

        for (uint64_t i = Qrange.fs, next_mega = 0; i < Qrange.fs + Qrange.s_sz && extracts.size() < vm["max_extract"].as<uint64_t>(); ++i)
        {
            std::deque<int64_t> extract;
            extract.push_front(i);
            while (next_mega++ < vm["max_mega"].as<uint64_t>() && next_megadot(extract))
            {
                if (extract.size() % 2 == 0)
                    pmegarange = &megadots_tail[megasources[extract[0]].first];
                else
                    pmegarange = &megadots_head[megasources[extract[0]].first];
                if (pmegarange->s_sz == 0)
                {
                    bool addflag = true;
                    if (vm["min_seg_num"].as<uint64_t>() > 0 && vm["max_seg_num"].as<uint64_t>() > 0 && (extract.size() < 2 * vm["min_seg_num"].as<uint64_t>() || extract.size() > 2 * vm["max_seg_num"].as<uint64_t>()))
                        addflag = false;
                    else
                    {
                        for (auto &extract_pre : extracts)
                            if (extract_dis(extract_pre, extract) < vm["diff_thres"].as<uint64_t>())
                            {
                                addflag = false;
                                break;
                            }
                    }
                    if (addflag)
                    {
                        extracts.push_back(extract);
                        if (extracts.size() == vm["max_extract"].as<uint64_t>())
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

    void source_until(uint64_t i)
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
        std::set<Dot *> visited;
        visited.insert(megasource.first);
        path.push_back(megasource.first);
        do
        {
            if (to_tail)
            {
                for (uint64_t i = 0; i < path[0]->s_sz; ++i)
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
    }

    bool next_dot(std::deque<Dot *> &path, std::set<Dot *> &visited, bool to_tail)
    {
        uint64_t is = 0;
        while (!path.empty())
        {
            for (uint64_t i = is; i < path[0]->s_sz; ++i)
            {
                Dot *pdot = sources[path[0]->fs + i];
                if (visited.count(pdot) || to_tail == (pdot->n < 0 || path[0]->n < 0))
                    continue;
                visited.insert(pdot);
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
        uint64_t is = 0;
        while (!extract.empty())
        {
            if (extract.size() % 2 == 0)
                pmegarange = &megadots_tail[megasources[extract[0]].first];
            else
                pmegarange = &megadots_head[megasources[extract[0]].first];
            for (uint64_t i = is; i < pmegarange->s_sz; ++i)
            {
                uint64_t idx = pmegarange->fs + i;
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
                is = extract[0] + 1 - pmegarange->fs;
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
        for (uint64_t i = 1; i < extract.size(); i += 2)
        {
            summin_score += edges[megasources[extract[i]].first->n]->min_score;
            sumdiffval += megasources[extract[i + 1]].first->val - megasources[extract[i]].first->val;
            if (new_node == edges[megasources[extract[i]].first->n]->head && sumdiffval < summin_score)
                return false;
        }
        return true;
    }

    uint64_t extract_dis(std::deque<int64_t> extract1, std::deque<int64_t> extract2)
    {
        uint64_t row = extract1.size() / 2 + 1, col = extract2.size() / 2 + 1;
        std::unique_ptr<uint64_t[]> PD(new uint64_t[row * col]), Ed(new uint64_t[col - 1]), Fd(new uint64_t[row - 1]);
        PD[0] = 0;
        for (uint64_t i = 1; i < row; ++i)
        {
            if (file2short.count(edges[megasources[extract1[2 * i - 1]].first->n]->name))
                Fd[i - 1] = megasources[extract1[2 * i - 1]].first->s - megasources[extract1[2 * i - 2]].first->s;
            else
                Fd[i - 1] = megasources[extract1[2 * i - 1]].first->lambda;
            PD[i] = PD[i - 1] + Fd[i - 1];
        }
        for (uint64_t j = 1; j < col; ++j)
        {
            if (file2short.count(edges[megasources[extract2[2 * j - 1]].first->n]->name))
                Ed[j - 1] = megasources[extract2[2 * j - 1]].first->s - megasources[extract2[2 * j - 2]].first->s;
            else
                Ed[j - 1] = megasources[extract2[2 * j - 1]].first->lambda;
            PD[j * row] = PD[(j - 1) * row] + Ed[j - 1];
        }
        for (uint64_t i = 1; i < row; ++i)
            for (uint64_t j = 1; j < col; ++j)
            {
                uint64_t F = PD[j * row + i - 1] + Fd[i - 1];
                uint64_t E = PD[(j - 1) * row + i] + Ed[j - 1];
                PD[j * row + i] = std::min(E, F);
                if (megasources[extract1[2 * i - 1]].first->n == megasources[extract2[2 * j - 1]].first->n)
                {
                    uint64_t segdiff;
                    if (file2short.count(edges[megasources[extract1[2 * i - 1]].first->n]->name))
                    {
                        int64_t inter = std::min(megasources[extract1[2 * i - 1]].first->s, megasources[extract2[2 * j - 1]].first->s) - std::max(megasources[extract1[2 * i - 2]].first->s, megasources[extract2[2 * j - 2]].first->s); // inter may be minus, so use signed type int64_t
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
        int64_t max_inter = 0; // max_inter need compare with int64_t inter, thereby being int64_t as well
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
        for (uint64_t i = 0; i < extracts.size(); ++i)
        {
            std::string align_query, align_mid, align_ref;
            for (uint64_t j = 0; j < extracts[i].size(); j += 2)
            {
                std::deque<Dot *> &path = megasources[extracts[i][j]].second;
                if (!align_query.empty())
                {
                    align_query.push_back(' ');
                    align_mid.push_back(' ');
                    align_ref.push_back(' ');
                }

                bool islocal = file2short.count(edges[path[0]->n]->name);
                std::pair<std::unique_ptr<uint8_t[]>, uint64_t> *ref = islocal ? (&file2short[edges[path[0]->n]->name]) : (&file2long[edges[path[0]->n]->name]);

                if (islocal)
                {
                    align_query.append(path[0]->s, '-');
                    align_mid.append(path[0]->s, ' ');
                    for (uint64_t i=0; i<path[0]->s; ++i)
                        align_ref.push_back(NUCTYPE2base[ref->first[i]]);
                }

                for (uint64_t i = 1; i < path.size(); ++i)
                {
                    bool E = path[i]->w > path[i - 1]->w, localF = islocal && path[i]->s > path[i - 1]->s, globalF = !islocal && path[i]->lambda > path[i - 1]->lambda;

                    if(!E && !localF && !globalF)
                        continue;

                    if (E)
                        align_query.push_back(O[path[i - 1]->w]);
                    else
                        align_query.push_back('-');
                    
                    if (localF)
                        align_ref.push_back(NUCTYPE2base[ref->first[path[i]->s - 1]]);
                    else if (globalF)
                        align_ref.push_back(NUCTYPE2base[ref->first[ref->second - 1 - path[i]->s]]);
                    else
                        align_ref.push_back('-');

                    if (std::tolower(align_ref.back()) == std::tolower(align_query.back()))
                        align_mid.push_back('|');
                    else
                        align_mid.push_back(' ');
                }

                if (islocal)
                {
                    uint64_t sz = ref->second - path.back()->s;
                    align_query.append(sz, '-');
                    align_mid.append(sz, ' ');
                    for (uint64_t i=path.back()->s; i<path.back()->s+sz; ++i)
                        align_ref.push_back(NUCTYPE2base[ref->first[i]]);
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
        for (uint64_t i = 0; i < extracts.size(); ++i)
        {
            fout_extract << Oname << '\t' << i << '\n';
            for (uint64_t j = 0; j < extracts[i].size(); j += 2)
            {
                Dot *pdot1 = megasources[extracts[i][j]].first;
                Dot *pdot2 = megasources[extracts[i][j + 1]].first;
                Edge *edge = edges[pdot2->n];
                if (file2rankvec.count(edge->name))
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
        uint8_t *longref = file2long[edges[pdot->n]->name].first.get();
        std::vector<std::pair<std::string, uint64_t>> &cumlen = file2cumlen[edges[pdot->n]->name];
        RankVec &rankvec = file2rankvec[edges[pdot->n]->name];
        int64_t sr1 = 0, sr2 = rankvec.bwt_sz;
 
        for (uint64_t s = 0, tmp = rankvec.bwt_sz - pdot->s + pdot->lambda - 2; s < pdot->lambda; ++s)
            rankvec.PreRange(sr1, sr2, longref[tmp - s]);
        
        std::ifstream &SAfin = file2SA[edges[pdot->n]->name];
        for (uint64_t sx = sr1; sx < sr2 && sx - sr1 < vm["max_range"].as<uint64_t>(); ++sx)
        {
            uint64_t s;
            SAfin.seekg(sx * sizeof(uint64_t));
            SAfin.read((char*)&s, sizeof(uint64_t));
            uint64_t l = 0, r = cumlen.size();
            while (r - l > 1)
            {
                uint64_t m = (l + r) / 2;
                if (cumlen[m].second > s)
                    r = m;
                else
                    l = m;
            }
            if (r < cumlen.size())
                s = cumlen[r].second - 1 - s;
            else
                s = rankvec.bwt_sz - 1 - s;
            seqranges[cumlen[l].first].emplace_back(s - pdot->lambda, s);
        }
        return seqranges;
    }
};

#endif