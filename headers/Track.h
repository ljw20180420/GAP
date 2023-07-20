#ifndef TRACK_H
#define TRACK_H

#include "Align.h"

struct Track : Graph
{
    struct MegaRange
    {
        SIZETYPE fs;
        SIZETYPE s_sz;
    };

    std::ofstream fout_align, fout_extract, fout_fail, fout_death;
    std::vector<Dot> dots;
    std::map<Dot *, MegaRange> megadots_head;
    std::map<Dot *, MegaRange> megadots_tail;
    std::deque<std::pair<Dot *, std::deque<Dot *>>> megasources;
    std::map<std::string, std::vector<std::pair<std::string, SIZETYPE>>> file2cumlen;

    Track(boost::program_options::variables_map &vm, std::map<std::string, std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE>> &file2short, std::map<std::string, RankVec> &file2rankvec)
    : Graph(vm, file2short, file2rankvec), fout_align("alg"), fout_extract("ext"), fout_fail("fail"), fout_death("death")
    {
        for (std::string file : vm["longs"].as<std::vector<std::string>>())
        {
            file = file.substr(0, file.find(','));
            if (file2cumlen.count(file))
                continue;

            std::string refname;
            SIZETYPE cumlen;
            for (std::ifstream fin(file+".cumlen"); fin >> refname >> cumlen;)
                file2cumlen[file].emplace_back(refname, cumlen);
        }
    }

    void ReadTrack(std::deque<std::string> mg_files)
    {
        std::vector<std::ifstream> fins;
        std::vector<SIZETYPE> fzs;
        std::vector<std::string> Onames;
        IDTYPE max_id = 0, id;
        for (std::string &mg_file : mg_files)
        {
            SIZETYPE fz = std::filesystem::file_size(mg_file);
            if (fz <= sizeof(IDTYPE))
                continue;
            fzs.push_back(fz);
            fins.emplace_back(mg_file, std::ios::binary);
            fins.back().seekg(-sizeof(IDTYPE), fins.back().end);
            fins.back().read((char *)&id, sizeof(IDTYPE));
            max_id = id > max_id ? id : max_id;
            fins.back().seekg(0, fins.back().beg);
            SIZETYPE Onsz;
            fins.back().read((char *)&Onsz, sizeof(Onsz));
            Onames.emplace_back(Onsz, 'x');
            fins.back().read((char *)Onames.back().data(), Onsz);
        }

        dots.resize(max_id + 1);
        std::ifstream finput(vm["input"].as<std::string>());
        for (std::string name, read; std::getline(std::getline(finput, name), read); )
        {
            for (SIZETYPE f = 0; f < fins.size(); ++f)
            {
                if (Onames[f] == name)
                {
                    std::deque<Dot *> sources;
                    for (IDTYPE id = 0, uid = 0; id <= uid; ++id)
                    {
                        fins[f].read((char *)&dots[id].n, sizeof(Dot::n));
                        fins[f].read((char *)&dots[id].s, sizeof(Dot::s));
                        fins[f].read((char *)&dots[id].w, sizeof(Dot::w));
                        fins[f].read((char *)&dots[id].val, sizeof(Dot::val));
                        fins[f].read((char *)&dots[id].s_sz, sizeof(Dot::s_sz));
                        fins[f].read((char *)&dots[id].lambda, sizeof(Dot::lambda));
                        dots[id].fs = sources.size();
                        for (SIZETYPE i = 0; i < dots[id].s_sz; ++i)
                        {
                            IDTYPE idd;
                            fins[f].read((char *)&idd, sizeof(idd));
                            sources.emplace_back(&dots[idd]);
                            uid = idd > uid ? idd : uid;
                        }
                    }

                    BackTrack(&dots[0], name, read, sources);

                    if (fzs[f]-fins[f].tellg()>sizeof(IDTYPE))
                    {
                        SIZETYPE Onsz;
                        fins[f].read((char *)&Onsz, sizeof(Onsz));
                        Onames[f].resize(Onsz);
                        fins[f].read((char *)Onames[f].data(), Onsz);
                    }

                    break;
                }
            }
        }
    }

    void BackTrack(Dot *pQ, std::string &Oname, std::string &O, std::deque<Dot *> &sources)
    {
        megasources.clear();
        megadots_head.clear();
        megadots_tail.clear();
        megasources.emplace_back(pQ, std::deque<Dot *>());

        for (SIZETYPE i = 0; i < megasources.size(); ++i)
            source_until(i, sources);
        MegaRange Qrange = megadots_tail[pQ];

        MegaRange *pmegarange;
        std::deque<std::deque<SIZETYPE>> extracts;

        for (SIZETYPE i = Qrange.fs, next_mega = 0; i < Qrange.fs + Qrange.s_sz && extracts.size() < vm["max_extract"].as<SIZETYPE>(); ++i)
        {
            std::deque<SIZETYPE> extract;
            extract.push_front(i);
            while (next_mega++ < vm["max_mega"].as<SIZETYPE>() && next_megadot(extract))
            {
                if (extract.size() % 2 == 0)
                    pmegarange = &megadots_tail[megasources[extract[0]].first];
                else
                    pmegarange = &megadots_head[megasources[extract[0]].first];
                if (pmegarange->s_sz == 0)
                {
                    bool addflag = true;
                    if (vm["min_seg_num"].as<SIZETYPE>() > 0 && vm["max_seg_num"].as<SIZETYPE>() > 0 && (extract.size() < 2 * vm["min_seg_num"].as<SIZETYPE>() || extract.size() > 2 * vm["max_seg_num"].as<SIZETYPE>()))
                        addflag = false;
                    else
                    {
                        for (std::deque<SIZETYPE> &extract_pre : extracts)
                            if (extract_dis(extract_pre, extract) < vm["diff_thres"].as<SIZETYPE>())
                            {
                                addflag = false;
                                break;
                            }
                    }
                    if (addflag)
                    {
                        extracts.push_back(extract);
                        if (extracts.size() == vm["max_extract"].as<SIZETYPE>())
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

    void source_until(SIZETYPE i, std::deque<Dot *> &sources)
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
                for (SIZETYPE i = 0; i < path[0]->s_sz; ++i)
                    if (sources[path[0]->fs + i]->n < 0)
                    {
                        megasources.emplace_back(path[0], path);
                        break;
                    }
            }
            else if (path[0]->n >= 0 && path[0] != megasource.first)
                megasources.emplace_back(path[0], std::deque<Dot *>());

        } while (next_dot(path, visited, to_tail, sources));
        insertpair.first->second.s_sz = megasources.size() - insertpair.first->second.fs;
    }

    bool next_dot(std::deque<Dot *> &path, std::set<Dot *> &visited, bool to_tail, std::deque<Dot *> &sources)
    {
        SIZETYPE is = 0;
        while (!path.empty())
        {
            for (SIZETYPE i = is; i < path[0]->s_sz; ++i)
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

    bool next_megadot(std::deque<SIZETYPE> &extract)
    {
        MegaRange *pmegarange;
        SIZETYPE is = 0;
        while (!extract.empty())
        {
            if (extract.size() % 2 == 0)
                pmegarange = &megadots_tail[megasources[extract[0]].first];
            else
                pmegarange = &megadots_head[megasources[extract[0]].first];
            for (SIZETYPE i = is; i < pmegarange->s_sz; ++i)
            {
                SIZETYPE idx = pmegarange->fs + i;
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

    bool legal_circle(Dot *pdot, std::deque<SIZETYPE> &extract)
    {
        if (extract.size() % 2 == 0)
            return true;
        SCORETYPE summin_score = edges[pdot->n]->min_score;
        SCORETYPE sumdiffval = megasources[extract[0]].first->val - pdot->val;
        Node *new_node = edges[pdot->n]->tail;
        if (new_node == edges[pdot->n]->head && sumdiffval < summin_score)
            return false;
        for (SIZETYPE i = 1; i < extract.size(); i += 2)
        {
            summin_score += edges[megasources[extract[i]].first->n]->min_score;
            sumdiffval += megasources[extract[i + 1]].first->val - megasources[extract[i]].first->val;
            if (new_node == edges[megasources[extract[i]].first->n]->head && sumdiffval < summin_score)
                return false;
        }
        return true;
    }

    SIZETYPE extract_dis(std::deque<SIZETYPE> extract1, std::deque<SIZETYPE> extract2)
    {
        SIZETYPE row = extract1.size() / 2 + 1, col = extract2.size() / 2 + 1;
        std::unique_ptr<SIZETYPE []> PD(new SIZETYPE[row * col]), Ed(new SIZETYPE[col - 1]), Fd(new SIZETYPE[row - 1]);
        PD[0] = 0;
        for (SIZETYPE i = 1; i < row; ++i)
        {
            if (file2short.count(edges[megasources[extract1[2 * i - 1]].first->n]->name))
                Fd[i - 1] = megasources[extract1[2 * i - 1]].first->s - megasources[extract1[2 * i - 2]].first->s;
            else
                Fd[i - 1] = megasources[extract1[2 * i - 1]].first->lambda;
            PD[i] = PD[i - 1] + Fd[i - 1];
        }
        for (SIZETYPE j = 1; j < col; ++j)
        {
            if (file2short.count(edges[megasources[extract2[2 * j - 1]].first->n]->name))
                Ed[j - 1] = megasources[extract2[2 * j - 1]].first->s - megasources[extract2[2 * j - 2]].first->s;
            else
                Ed[j - 1] = megasources[extract2[2 * j - 1]].first->lambda;
            PD[j * row] = PD[(j - 1) * row] + Ed[j - 1];
        }
        for (SIZETYPE i = 1; i < row; ++i)
            for (SIZETYPE j = 1; j < col; ++j)
            {
                SIZETYPE F = PD[j * row + i - 1] + Fd[i - 1];
                SIZETYPE E = PD[(j - 1) * row + i] + Ed[j - 1];
                PD[j * row + i] = std::min(E, F);
                if (megasources[extract1[2 * i - 1]].first->n == megasources[extract2[2 * j - 1]].first->n)
                {
                    SIZETYPE segdiff;
                    if (file2short.count(edges[megasources[extract1[2 * i - 1]].first->n]->name))
                    {
                        SIZETYPE right_min = std::min(megasources[extract1[2 * i - 1]].first->s, megasources[extract2[2 * j - 1]].first->s);
                        SIZETYPE left_max = std::max(megasources[extract1[2 * i - 2]].first->s, megasources[extract2[2 * j - 2]].first->s);
                        SIZETYPE inter = right_min > left_max ? right_min - left_max : 0;
                        segdiff = Fd[i - 1] + Ed[j - 1] - 2 * inter;
                    }
                    else
                    {
                        std::map<std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> seqranges1 = get_seqranges(megasources[extract1[2 * i - 1]].first), seqranges2 = get_seqranges(megasources[extract2[2 * j - 1]].first);
                        segdiff = global_dis(seqranges1, seqranges2);
                    }
                    PD[j * row + i] = std::min(PD[j * row + i], PD[(j - 1) * row + i - 1] + segdiff);
                }
            }
        return PD[row * col - 1];
    }

    SIZETYPE global_dis(std::map<std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> &seqranges1, std::map<std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> &seqranges2)
    {
        SIZETYPE max_inter = 0;
        for (auto it1 = seqranges1.begin(), it2 = seqranges2.begin(); it1 != seqranges1.end() && it2 != seqranges2.end();)
        {
            int cmp = it1->first.compare(it2->first);
            if (cmp < 0)
                ++it1;
            else if (cmp > 0)
                ++it2;
            else
            {
                for (std::pair<SIZETYPE, SIZETYPE> &range1 : it1->second)
                    for (std::pair<SIZETYPE, SIZETYPE> &range2 : it2->second)
                    {
                        SIZETYPE right_min = std::min(range1.second, range2.second), left_max = std::max(range1.first, range2.first);
                        max_inter = right_min > left_max ? std::max(max_inter, right_min - left_max) : max_inter;
                    }
                ++it1;
                ++it2;
            }
        }
        std::pair<SIZETYPE, SIZETYPE> &range1 = (seqranges1.begin())->second.front(), &range2 = (seqranges2.begin())->second.front();
        return range1.second - range1.first + range2.second - range2.first - 2 * max_inter;
    }

    void extracts_to_aligns(std::deque<std::deque<SIZETYPE>> &extracts, std::string &Oname, std::string &O)
    {
        for (SIZETYPE i = 0; i < extracts.size(); ++i)
        {
            std::string align_query, align_mid, align_ref;
            for (SIZETYPE j = 0; j < extracts[i].size(); j += 2)
            {
                std::deque<Dot *> &path = megasources[extracts[i][j]].second;
                if (!align_query.empty())
                {
                    align_query.push_back(' ');
                    align_mid.push_back(' ');
                    align_ref.push_back(' ');
                }

                bool islocal = file2short.count(edges[path[0]->n]->name);
                NUCTYPE *ref;
                SIZETYPE ref_sz;
                NUCTYPE revref[path.back()->lambda];
                if (islocal)
                {
                    std::pair<std::unique_ptr<NUCTYPE []>, SHORTSIZE> &pair = file2short[edges[path[0]->n]->name];
                    ref = pair.first.get();
                    ref_sz = pair.second;
                }
                else
                {
                    std::ifstream &REVREFfin = file2long[edges[path[0]->n]->name];
                    REVREFfin.seekg(file2rankvec[edges[path[0]->n]->name].bwt_sz - 1 - path.back()->s);
                    REVREFfin.read((char*)revref, path.back()->lambda);
                }

                if (islocal)
                {
                    align_query.append(path[0]->s, '-');
                    align_mid.append(path[0]->s, ' ');
                    for (SIZETYPE i=0; i<path[0]->s; ++i)
                        align_ref.push_back(NUCTYPE2base[ref[i]]);
                }

                for (SIZETYPE i = 1; i < path.size(); ++i)
                {
                    bool E = path[i]->w > path[i - 1]->w, localF = islocal && path[i]->s > path[i - 1]->s, globalF = !islocal && path[i]->lambda > path[i - 1]->lambda;

                    if(!E && !localF && !globalF)
                        continue;

                    if (E)
                        align_query.push_back(O[path[i - 1]->w]);
                    else
                        align_query.push_back('-');
                    
                    if (localF)
                        align_ref.push_back(NUCTYPE2base[ref[path[i]->s - 1]]);
                    else if (globalF)
                        align_ref.push_back(NUCTYPE2base[revref[path.back()->lambda - path[i]->lambda]]);
                    else
                        align_ref.push_back('-');

                    if (std::tolower(align_ref.back()) == std::tolower(align_query.back()))
                        align_mid.push_back('|');
                    else
                        align_mid.push_back(' ');
                }

                if (islocal)
                {
                    SIZETYPE sz = ref_sz - path.back()->s;
                    align_query.append(sz, '-');
                    align_mid.append(sz, ' ');
                    for (SIZETYPE i=path.back()->s; i<path.back()->s+sz; ++i)
                        align_ref.push_back(NUCTYPE2base[ref[i]]);
                }
            }

            fout_align << Oname << '\t' << i << '\n'
                       << align_query << '\n'
                       << align_mid << '\n'
                       << align_ref << '\n';
        }
    }

    void output_extracts(std::deque<std::deque<SIZETYPE>> &extracts, std::string &Oname)
    {
        for (SIZETYPE i = 0; i < extracts.size(); ++i)
        {
            fout_extract << Oname << '\t' << i << '\n';
            for (SIZETYPE j = 0; j < extracts[i].size(); j += 2)
            {
                Dot *pdot1 = megasources[extracts[i][j]].first;
                Dot *pdot2 = megasources[extracts[i][j + 1]].first;
                Edge *edge = edges[pdot2->n];
                if (file2rankvec.count(edge->name))
                {
                    std::map<std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> seqranges = get_seqranges(pdot2);
                    fout_extract << edge->name << ':';
                    for (std::pair<const std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> &seqrange : seqranges)
                        for (std::pair<SIZETYPE, SIZETYPE> &range : seqrange.second)
                            fout_extract << seqrange.first << '\t' << range.first << '\t';
                    fout_extract << pdot1->w << '\t' << pdot1->val << '\n';
                    fout_extract << edge->name << ':';
                    for (std::pair<const std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> &seqrange : seqranges)
                        for (std::pair<SIZETYPE, SIZETYPE> &range : seqrange.second)
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

    std::map<std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> get_seqranges(Dot *pdot)
    {
        std::map<std::string, std::deque<std::pair<SIZETYPE, SIZETYPE>>> seqranges;
        std::vector<std::pair<std::string, SIZETYPE>> &cumlen = file2cumlen[edges[pdot->n]->name];
        RankVec &rankvec = file2rankvec[edges[pdot->n]->name];
        std::ifstream &REVREFfin = file2long[edges[pdot->n]->name];
        NUCTYPE revref[pdot->lambda];
        REVREFfin.seekg(rankvec.bwt_sz - pdot->s - 1);
        REVREFfin.read((char*)revref, pdot->lambda);

        SIZETYPE sr1 = 0, sr2 = rankvec.bwt_sz;
        for (SIZETYPE s = 1; s <= pdot->lambda; ++s)
        {
            NUCTYPE c = revref[pdot->lambda - s];
            sr1 = rankvec.C[c] + rankvec.rank(c,sr1);
            sr2 = rankvec.C[c] + rankvec.rank(c,sr2);
        }
            
        
        std::ifstream &SAfin = file2SA[edges[pdot->n]->name];
        for (SIZETYPE sx = sr1; sx < sr2 && sx - sr1 < vm["max_range"].as<SIZETYPE>(); ++sx)
        {
            SIZETYPE s;
            SAfin.seekg(sx * sizeof(SIZETYPE));
            SAfin.read((char*)&s, sizeof(SIZETYPE));
            SIZETYPE l = 0, r = cumlen.size();
            while (r - l > 1)
            {
                SIZETYPE m = (l + r) / 2;
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