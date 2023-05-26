#ifndef TRACK_H
#define TRACK_H

#include "Align.h"

struct Track : Graph
{
    MonoDeque<Dot *> sources;

    Track(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel) : Graph(argc, argv, file2seq, file2browheel)
    {
    }

    void ReadTrack(std::deque<std::string> mg_files, std::deque<std::string> read_files, std::string run_name, int64_t max_seq, int64_t max_track)
    {
        std::vector<std::ifstream> fins;
        for (auto &file : mg_files)
            fins.emplace_back(file, std::ios::binary);
        std::ofstream fout_track(run_name+".trk");
        std::ofstream fout_align(run_name+".alg");
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
                fins[i].read((char*)Onames[i].data(), sizeof(char) * Onsz);
            }
        }
        std::vector<Dot> dots(max_id + 1);
        int64_t seq_num = 0;
        std::string Oname, O;
        auto iter = read_files.begin();
        std::ifstream fin_seq(*iter);
        while (seq_num < max_seq && iter != read_files.end())
        {
            if (std::getline(std::getline(fin_seq, Oname), O))
            {
                ++seq_num;
                if (Oname[0]=='@')
                    fin_seq.ignore(std::numeric_limits<std::streamsize>::max(),'\n').ignore(std::numeric_limits<std::streamsize>::max(),'\n');
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

                        std::deque<Dot *> path;
                        path.push_front(&dots[0]);
                        BackTrack(fout_track, fout_align, path, max_track, Oname, O);

                        sources.clear();
                        if (fzs[f] - fins[f].tellg() > sizeof(int64_t))
                        {
                            int Onsz;
                            fins[f].read((char *)&Onsz, sizeof(Onsz));
                            Onames[f].resize(Onsz);
                            fins[f].read((char*)Onames[f].data(), sizeof(char) * Onsz);
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

    void BackTrack(std::ofstream &fout_track, std::ofstream &fout_align, std::deque<Dot *> &path, int64_t max_track, std::string &Oname, std::string &O)
    {
        int64_t track_num = 0;
        while (!path.empty() && track_num < max_track)
        {
            if (path[0]->s_sz == 0)
            {
                int node_n = path[1]->n;
                fout_track << Oname << '\t' << ++track_num << '\t' << path.size() << '\n' << node_n << '\t' << path[0]->s << '\t' << path[0]->w << '\t' << path[0]->val << '\n';
                for (int i = 1; i<path.size(); ++i)
                {
                    if (path[i]->n >= 0)
                    {
                        if (edges[path[i]->n]->global)
                        {
                            int n = path[i]->n;
                            int is = i;
                            do
                            {
                                ++i;
                            } while (i < path.size() && path[i]->n == n);
                            int lambda = path[i-1]->lambda;
                            int64_t s = path[i-1]->s, ls;
                            std::string name;
                            do
                            {
                                std::tie(name, ls) = static_cast<EdgeGlobal *>(edges[path[is]->n])->browheel.get_axis(s + path[is]->lambda - lambda);
                                fout_track << path[is]->n << '\t' << name << '\t' << ls << '\t' << path[is]->w << '\t' << path[is]->val << '\n';
                                ++is;
                            } while (is != i);
                            --i;
                        }
                        else
                            fout_track << path[i]->n << '\t' << path[i]->s << '\t' << path[i]->w << '\t' << path[i]->val << '\n';
                    }
                    else
                    {
                        int n = path[i-1]->n;
                        if (n >= 0)
                            node_n = Dot::nidx_trans(edges[n]->head->n);
                        fout_track << node_n << '\t' << path[i]->s << '\t' << path[i]->w << '\t' << path[i]->val << '\n';
                    }
                }

                std::string align_query, align_mid, align_ref;
                for (int i = 1; i < path.size() - 1; ++i)
                {
                    if (path[i]->w > path[i - 1]->w && (path[i]->s > path[i - 1]->s || (path[i]->n == path[i - 1]->n && path[i]->n >= 0 && edges[path[i]->n]->global && path[i]->lambda > path[i - 1]->lambda)))
                    {
                        char c;
                        if (!edges[path[i]->n]->global)
                            c = static_cast<EdgeLocal *>(edges[path[i]->n])->nameseq.seq[path[i]->s - 1];
                        else
                        {
                            auto &browheel = static_cast<EdgeGlobal *>(edges[path[i]->n])->browheel;
                            c = BroWheel::int2base[browheel.sequence(browheel.sequence.size() - 1 - path[i]->s)];
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
                                int sz = (path[i]->n < 0 ? local->nameseq.seq.size() : path[i]->s) - (path[i - 1]->n < 0 ? 0 : path[i - 1]->s);
                                align_query.append(sz, '-');
                                align_mid.append(sz, ' ');
                                align_ref.append(local->nameseq.seq, path[i - 1]->s, sz);
                            }
                            else
                            {
                                auto &browheel = static_cast<EdgeGlobal *>(edges[n])->browheel;
                                align_query.push_back('-');
                                align_mid.push_back(' ');
                                align_ref.push_back(BroWheel::int2base[browheel.sequence(browheel.sequence.size() - 1 - path[i]->s)]);
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
                fout_align << Oname << '\t' << track_num << '\n' << align_query << '\n' << align_mid << '\n' << align_ref << '\n';

                Dot *dot_parent = *std::next(path.begin());
                while (path.front() == sources[dot_parent->fs + dot_parent->s_sz - 1])
                {
                    path.pop_front();
                    if (path.size() == 1)
                    {
                        path.pop_front();
                        break;
                    }
                    dot_parent = *std::next(path.begin());
                }
                if (!path.empty())
                    for (int i = 0; i < dot_parent->s_sz; ++i)
                        if (sources[dot_parent->fs + i] == path.front())
                        {
                            path.pop_front();
                            path.push_front(sources[path.front()->fs + i + 1]);
                            break;
                        }
            }
            else
                path.push_front(sources[path.front()->fs]);
        }
    }

    void extract(std::string run_name, int64_t max_extract)
    {
        std::ifstream fin(run_name + ".trk");
        std::ofstream fout(run_name + ".ext");
        std::string Oname, line;
        int64_t track_num;
        int64_t track_length;
        while (true)
        {
            do
            {
                std::getline(fin, line);
                std::istringstream iss(line);
                std::getline(iss, Oname, '\t');
                iss >> track_num >> track_length;
                if (!fin.good() || track_num <= max_extract)
                    break;
                for (int64_t i = 0; i < track_length; ++i)
                    std::getline(fin, line);
            } while (true);
            if (!fin.good())
                break;
            fout << Oname << '\t' << track_num << '\n';
            if (track_length == 0)
                continue;
            else
            {

                int np, wp, n, w;
                int64_t sp, s;
                double valp, val;
                std::string namep, name;
                fin >> n >> s >> w >> val >> std::ws;
                for (int64_t i = 1; i < track_length; ++i)
                {
                    np = n;
                    wp = w;
                    sp = s;
                    valp = val;
                    namep.swap(name);
                    fin >> n;
                    if (n >= 0 && edges[n]->global)
                        fin >> name;
                    fin >> s >> w >> val >> std::ws;
                    if (np != n)
                    {
                        if (np < 0)
                        {
                            fout << edges[n]->name;
                            if (edges[n]->global)
                                fout << ':' << name;
                            fout << '\t' << s << '\t' << w << '\t' << val << '\n';
                        }
                        else
                        {
                            fout << edges[np]->name;
                            if (edges[np]->global)
                                fout << ':' << namep;
                            fout << '\t' << sp << '\t' << wp << '\t' << valp << '\n';
                        }
                    }
                }
            }
        }
    }
};

#endif