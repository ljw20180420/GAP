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

    struct MegaDot
    {
        Dot *pdot;
        int fs;
        int s_sz;
    };

    MonoDeque<Dot *> sources;
    MonoDeque<std::pair<MegaDot *, std::deque<Dot *>>> mega_sources;

    Track(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel) : Graph(argc, argv, file2seq, file2browheel)
    {
    }

    void ReadTrack(std::deque<std::string> mg_files, std::deque<std::string> read_files, std::string run_name, int64_t max_seq, int64_t max_track, int64_t max_extract, int round)
    {
        std::vector<std::ifstream> fins;
        for (auto &file : mg_files)
            fins.emplace_back(file, std::ios::binary);
        std::ofstream fout_track(run_name + ".trk" + std::to_string(round)), fout_align(run_name + ".alg" + std::to_string(round)), fout_extract(run_name + ".ext" + std::to_string(round)), fout_fail(run_name + ".fail" + std::to_string(round));
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
        std::vector<Dot> dots(max_id + 1);
        int64_t seq_num = 0;
        std::string Oname, O, plr, quality;
        auto iter = read_files.begin();
        std::ifstream fin_seq(*iter);
        while (seq_num < max_seq && iter != read_files.end())
        {
            if (std::getline(std::getline(fin_seq, Oname), O))
            {
                ++seq_num;
                if (Oname[0] == '@')
                    std::getline(std::getline(fin_seq, plr), quality);
                for (int64_t f = 0; f < fins.size(); ++f)
                {
                    if (Onames[f] == Oname)
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

                        BackTrack(fout_track, fout_align, fout_extract, fout_fail, &dots[0], max_track, max_extract, Oname, O, plr, quality);

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
            else
            {
                fin_seq.close();
                fin_seq.open(*(++iter));
            }
        }
    }

    void BackTrack(std::ofstream &fout_track, std::ofstream &fout_align, std::ofstream &fout_extract, std::ofstream &fout_fail, Dot *pQ, int64_t max_track, int64_t max_extract, std::string &Oname, std::string &O, std::string &plr, std::string &quality)
    {
        std::deque<std::deque<Peak>> extracts;
        std::deque<Peak> extract;
        Peak peak;
        int64_t track_num = 0;
        std::deque<Dot *> path{pQ};
        while (!path.empty() && track_num < max_track)
        {
            if (next_dot(path, peak, extract, extracts) && path[0]->s_sz == 0)
            {
                extracts.push_back(extract);
                int n = path[1]->n;
                fout_track << Oname << '\t' << ++track_num << '\t' << path.size() << '\n'
                           << n << '\t' << path[0]->s << '\t' << path[0]->w << '\t' << path[0]->val << '\n';
                for (int i = 1, e = 0; i < path.size(); ++i)
                {
                    if ((path[i]->n >= 0) == (path[i - 1] < 0))
                    {
                        n = path[i]->n;
                        if (path[i]->n < 0 && path[i - 1]->n >= 0)
                            ++e;
                    }
                    if (path[i]->n >= 0 && edges[path[i]->n]->global)
                        fout_track << extract[e].n << '\t' << extract[e].seqname << '\t' << extract[e].sranges[0].first + path[i]->lambda << '\t' << path[i]->w << '\t' << path[i]->val << '\n';
                    else
                        fout_track << n << '\t' << path[i]->s << '\t' << path[i]->w << '\t' << path[i]->val << '\n';
                }

                std::string align_query, align_mid, align_ref;
                for (int i = 1; i < path.size() - 1; ++i)
                {
                    if (path[i]->w > path[i - 1]->w && (path[i]->s > path[i - 1]->s || (path[i]->n == path[i - 1]->n && path[i]->n >= 0 && edges[path[i]->n]->global && path[i]->lambda > path[i - 1]->lambda)))
                    {
                        char c;
                        if (!edges[path[i]->n]->global)
                            c = static_cast<EdgeLocal *>(edges[path[i]->n])->pnameseq->seq[path[i]->s - 1];
                        else
                        {
                            BroWheel *pbrowheel = static_cast<EdgeGlobal *>(edges[path[i]->n])->pbrowheel;
                            c = BroWheel::int2base[pbrowheel->sequence(pbrowheel->sequence.size() - 1 - path[i]->s)];
                        }
                        align_query.push_back(O[path[i]->w - 1]);
                        align_mid.push_back('|');
                        align_ref.push_back(c);
                    }
                    else
                    {
                        int sz = path[i]->w - path[i - 1]->w;
                        align_query.append(O, path[i - 1]->w, sz);
                        align_mid.append(sz, ' ');
                        align_ref.append(sz, '-');

                        if (path[i - 1]->n < 0 && path[i]->n >= 0)
                        {
                            align_query.push_back(' ');
                            align_mid.push_back(' ');
                            align_ref.push_back(' ');
                        }

                        int n = std::max(path[i - 1]->n, path[i]->n);
                        if (n >= 0 && (std::min(path[i - 1]->n, path[i]->n) >= 0 || !edges[n]->global))
                        {
                            if (!edges[n]->global)
                            {
                                EdgeLocal *local = static_cast<EdgeLocal *>(edges[n]);
                                int sz = (path[i]->n < 0 ? local->pnameseq->seq.size() : path[i]->s) - (path[i - 1]->n < 0 ? 0 : path[i - 1]->s);
                                align_query.append(sz, '-');
                                align_mid.append(sz, ' ');
                                align_ref.append(local->pnameseq->seq, path[i - 1]->s, sz);
                            }
                            else if (path[i]->s != path[i - 1]->s)
                            {
                                BroWheel *pbrowheel = static_cast<EdgeGlobal *>(edges[n])->pbrowheel;
                                align_query.push_back('-');
                                align_mid.push_back(' ');
                                align_ref.push_back(BroWheel::int2base[pbrowheel->sequence(pbrowheel->sequence.size() - 1 - path[i]->s)]);
                            }
                        }
                        if (path[i - 1]->n >= 0 && path[i]->n < 0)
                        {
                            align_query.push_back(' ');
                            align_mid.push_back(' ');
                            align_ref.push_back(' ');
                        }
                    }
                }
                fout_align << Oname << '\t' << track_num << '\n'
                           << align_query << '\n'
                           << align_mid << '\n'
                           << align_ref << '\n';
            }
        }
        for (int i = 0, j = 0; j < max_extract && i < extracts.size(); ++i)
        {
            std::vector<std::string> edgeseqnames(extracts[i].size());
            for (int k = 0; k < extracts[i].size(); ++k)
            {
                edgeseqnames[k] = edges[extracts[i][k].n]->name;
                if (!extracts[i][k].seqname.empty())
                    edgeseqnames[k] += ":" + extracts[i][k].seqname;
            }
            std::vector<int64_t> idxs(extracts[i].size(), 0);
            int kk;
            do
            {
                fout_extract << Oname << '\t' << i + 1 << '\t' << ++j << '\n';
                for (int k = 0; k < idxs.size(); ++k)
                {
                    Peak &peak_now = extracts[i][k];
                    auto &srange = peak_now.sranges[idxs[k]];
                    fout_extract << edgeseqnames[k] << '\t' << srange.first << '\t' << peak_now.wrange.first << '\t' << peak_now.valrange.first << '\n'
                                 << edgeseqnames[k] << '\t' << srange.second << '\t' << peak_now.wrange.second << '\t' << peak_now.valrange.second << '\n';
                }
                for (kk = 0; kk < idxs.size(); ++kk)
                {
                    ++idxs[kk];
                    if (idxs[kk] == extracts[i][kk].sranges.size())
                        idxs[kk] = 0;
                    else
                        break;
                }
            } while (j < max_extract && kk < idxs.size());
        }
        if (track_num == 0)
            fout_fail << Oname << '\n'
                      << O << '\n'
                      << plr << '\n'
                      << quality << '\n';
    }

    bool next_dot(std::deque<Dot *> &path, Peak &peak, std::deque<Peak> &extract, std::deque<std::deque<Peak>> &extracts)
    {
        int is = 0;
        while (!path.empty())
        {
            for (int i = is; i < path[0]->s_sz; ++i)
            {
                Dot *dotp = sources[path[0]->fs + i];
                if (dotp->n >= 0 && path[0]->n < 0)
                {
                    peak.n = dotp->n;
                    if (edges[dotp->n]->global)
                    {
                        BroWheel *pbrowheel = static_cast<EdgeGlobal *>(edges[dotp->n])->pbrowheel;
                        int64_t sr1 = 0, sr2 = pbrowheel->sequence.size() - 1;
                        for (int64_t s = pbrowheel->sequence.size() - dotp->s + dotp->lambda - 2; s >= pbrowheel->sequence.size() - 1 - dotp->s; --s)
                            pbrowheel->PreRange(sr1, sr2, pbrowheel->sequence(s));
                        peak.sranges.clear();
                        for (int64_t sx = sr1; sx <= sr2; ++sx)
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
                            if (sx == sr1)
                                peak.seqname = pbrowheel->name_cumlen[l].first;
                            if (r < pbrowheel->name_cumlen.size())
                                s = pbrowheel->name_cumlen[r].second - 1 - s;
                            else
                                s = pbrowheel->sequence.size() - 1 - s;
                            peak.sranges.emplace_back(s - dotp->lambda, s);
                        }
                        filter_sranges(peak, extract, extracts, edges[dotp->n]->diffseg);
                        if (peak.sranges.empty())
                            continue;
                        if (!legal_circle(peak, extract))
                            continue;
                    }
                    else
                    {
                        peak.seqname.clear();
                        peak.sranges.resize(1);
                        peak.sranges[0].second = dotp->s;
                    }
                    peak.wrange.second = dotp->w;
                    peak.valrange.second = dotp->val;
                }
                if (dotp->n < 0 && path[0]->n >= 0)
                {
                    if (!edges[path[0]->n]->global)
                    {
                        peak.sranges[0].first = path[0]->s;
                        filter_sranges(peak, extract, extracts, edges[path[0]->n]->diffseg);
                        if (peak.sranges.empty())
                            continue;
                        if (!legal_circle(peak, extract))
                            continue;
                    }
                    peak.wrange.first = path[0]->w;
                    peak.valrange.first = path[0]->val;
                    extract.push_front(peak);
                }
                path.push_front(dotp);
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
                if (path[0]->n < 0 && path[1]->n >= 0)
                {
                    peak = extract[0];
                    extract.pop_front();
                }
            }
            path.pop_front();
        }

        return false;
    }

    void filter_sranges(Peak &peak, std::deque<Peak> &extract, std::deque<std::deque<Peak>> &extracts, int64_t diffseg)
    {
        std::deque<std::pair<int64_t, int64_t>> sranges;
        sranges.swap(peak.sranges);
        for (auto &srange : sranges)
        {
            bool dup = false;
            for (auto &srange_pre : peak.sranges)
            {
                int64_t dis = std::abs(srange.first - srange_pre.first) + std::abs(srange.second - srange_pre.second);
                if (dis < diffseg)
                {
                    dup = true;
                    break;
                }
            }
            if (dup)
                continue;
            for (std::deque<Peak> &extract_pre : extracts)
            {
                if (extract_pre.size() <= extract.size())
                    continue;
                Peak &peak_pre = extract_pre[extract.size()];
                if (peak_pre.n != peak.n || peak_pre.seqname != peak.seqname)
                    continue;
                for (auto &srange_pre : peak_pre.sranges)
                {
                    int64_t dis = std::abs(srange.first - srange_pre.first) + std::abs(srange.second - srange_pre.second);
                    if (dis < diffseg)
                    {
                        dup = true;
                        break;
                    }
                }
                if (dup)
                    break;
            }
            if (dup)
                continue;
            peak.sranges.push_back(srange);
        }
    }

    bool legal_circle(Peak &peak, std::deque<Peak> &extract)
    {
        double summinScore = edges[peak.n]->minScore;
        double sumdiffval = peak.valrange.second - peak.valrange.first;
        Node *new_node = edges[peak.n]->tail;
        if (new_node == edges[peak.n]->head && sumdiffval < summinScore)
            return false;
        for (Peak &peak_pre : extract)
        {
            summinScore += edges[peak_pre.n]->minScore;
            sumdiffval += peak_pre.valrange.second - peak_pre.valrange.first;
            if (new_node == edges[peak_pre.n]->head && sumdiffval < summinScore)
                return false;
        }
        return true;
    }

    void source_until(Dot *pdot, std::string mode)
    {
        std::deque<Dot *> path;
        std::queue<Dot *> visited;
        pdot->visit = true;
        
    }
};

#endif