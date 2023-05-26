#ifndef TRACK_H
#define TRACK_H

#include "Align.h"

struct Track : Memory, Graph
{
    std::vector<EdgeLocal *> locals;
    std::vector<EdgeGlobal *> globals;
    std::set<std::string> node_names;

    Track(int argc, char **argv, std::map<std::string, NameSeq> &file2seq, std::map<std::string, BroWheel> &file2browheel, size_t chunk_sz_) : Graph(argc, argv, file2seq, file2browheel)
    {
        Initial(chunk_sz_);
        for (auto &node : nodes)
            node_names.emplace(node.name);
        int max_n = 0;
        for (auto &edge : local_crosses)
        {
            if (edge.n > max_n)
                max_n = edge.n;
            if (size_t(max_n) >= locals.size())
                locals.resize(max_n + 1, NULL);
            locals[edge.n] = &edge;
        }
        for (auto &edge : local_circuits)
        {
            if (edge.n > max_n)
                max_n = edge.n;
            if (size_t(max_n) >= locals.size())
                locals.resize(max_n + 1, NULL);
            locals[edge.n] = &edge;
        }
        for (auto &edge : global_crosses)
        {
            if (edge.n > max_n)
                max_n = edge.n;
            if (size_t(max_n) >= globals.size())
                globals.resize(max_n + 1, NULL);
            globals[edge.n] = &edge;
        }
        for (auto &edge : global_circuits)
        {
            if (edge.n > max_n)
                max_n = edge.n;
            if (size_t(max_n) >= globals.size())
                globals.resize(max_n + 1, NULL);
            globals[edge.n] = &edge;
        }
        locals.resize(max_n + 1, NULL);
        globals.resize(max_n + 1, NULL);
    }

    void ReadTrack(std::deque<std::string> mg_files, std::deque<std::string> read_files, std::string run_name, size_t max_seq, size_t max_track)
    {
        std::vector<std::ifstream> fins;
        for (auto &file : mg_files)
            fins.emplace_back(file, std::ios::binary);
        std::ofstream fout_track(run_name+".trk");
        std::ofstream fout_align(run_name+".alg");
        size_t max_id = 0;
        std::vector<std::string> Onames(fins.size());
        std::vector<size_t> fzs(fins.size());
        for (size_t i = 0; i < fins.size(); ++i)
        {
            size_t id;
            fins[i].seekg(0, fins[i].end);
            fzs[i] = fins[i].tellg();
            if (fzs[i] > sizeof(size_t))
            {
                fins[i].seekg(-sizeof(size_t), fins[i].end);
                fins[i].read((char *)&id, sizeof(size_t));
                max_id = id > max_id ? id : max_id;
                fins[i].seekg(0, fins[i].beg);
                int Onsz;
                fins[i].read((char *)&Onsz, sizeof(Onsz));
                Onames[i].resize(Onsz);
                fins[i].read(Onames[i].data(), sizeof(char) * Onsz);
            }
        }
        std::vector<Dot> dots(max_id + 1);
        size_t seq_num = 0;
        std::string Oname, O;
        auto iter = read_files.begin();
        std::ifstream fin_seq(*iter);
        while (seq_num < max_seq && iter != read_files.end())
        {
            if (std::getline(std::getline(fin_seq, Oname), O))
            {
                ++seq_num;
                if (Oname[0]=='>')
                    fin_seq.ignore(std::numeric_limits<size_t>::max(),'\n').ignore(std::numeric_limits<size_t>::max(),'\n');
                for (size_t f = 0; f < fins.size(); ++f)
                {
                    if (Onames[f] == Oname)
                    {
                        ssize_t id = -1;
                        size_t uid = 0;
                        do
                        {
                            ++id;
                            fins[f].read((char *)&dots[id].n, sizeof(Dot::n));
                            fins[f].read((char *)&dots[id].s, sizeof(Dot::s));
                            fins[f].read((char *)&dots[id].w, sizeof(Dot::w));
                            fins[f].read((char *)&dots[id].val, sizeof(Dot::val));
                            fins[f].read((char *)&dots[id].s_sz, sizeof(Dot::s_sz));
                            fins[f].read((char *)&dots[id].lambda, sizeof(Dot::lambda));
                            dots[id].sources = heap_alloc<Dot *>(dots[id].s_sz);
                            for (int i = 0; i < dots[id].s_sz; ++i)
                            {
                                size_t idd;
                                fins[f].read((char *)&idd, sizeof(idd));
                                dots[id].sources[i] = &dots[idd];
                                uid = idd > uid ? idd : uid;
                            }
                        } while (size_t(id) < uid);

                        std::deque<Dot *> path;
                        path.push_front(&dots[0]);
                        BackTrack(fout_track, fout_align, path, max_track, Oname, O);

                        clear();
                        if (fzs[f] - fins[f].tellg() > sizeof(size_t))
                        {
                            int Onsz;
                            fins[f].read((char *)&Onsz, sizeof(Onsz));
                            Onames[f].resize(Onsz);
                            fins[f].read(Onames[f].data(), sizeof(char) * Onsz);
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

    void BackTrack(std::ofstream &fout_track, std::ofstream &fout_align, std::deque<Dot *> &path, size_t max_track, std::string &Oname, std::string &O)
    {
        size_t track_num = 0;
        while (!path.empty() && track_num < max_track)
        {
            if (path.front()->s_sz == 0)
            {
                fout_track << Oname << '\t' << ++track_num << '\t' << path.size() << '\n';
                std::string node_name = nodes[-(*std::next(path.begin()))->n - 4].name;
                for (auto iter = path.begin(); iter != path.end(); ++iter)
                {
                    if ((*iter)->n >= 0)
                    {
                        if (globals[(*iter)->n])
                        {
                            int n = (*iter)->n;
                            auto iter_s = iter;
                            do
                            {
                                ++iter;
                            } while (iter != path.end() && (*iter)->n == n);
                            int lambda = (*std::prev(iter))->lambda;
                            size_t s = (*std::prev(iter))->s, ls;
                            std::string name;
                            do
                            {
                                std::tie(name, ls) = globals[(*iter_s)->n]->browheel.get_axis(s + (*iter_s)->lambda - lambda);
                                fout_track << globals[(*iter_s)->n]->name << ':' << name << '\t' << ls << '\t' << (*iter_s)->w << '\t' << (*iter_s)->val << '\n';
                                ++iter_s;
                            } while (iter_s != iter);
                            --iter;
                        }
                        else
                        {
                            fout_track << locals[(*iter)->n]->name << ':' << locals[(*iter)->n]->nameseq.name << '\t' << (*iter)->s << '\t' << (*iter)->w << '\t' << (*iter)->val << '\n';
                        }
                    }
                    else
                    {
                        if (iter != path.begin())
                        {
                            int n = (*std::prev(iter))->n;
                            if (n >= 0)
                                node_name = locals[n] ? locals[n]->head->name : globals[n]->head->name;
                        }
                        fout_track << node_name << '\t' << (*iter)->s << '\t' << (*iter)->w << '\t' << (*iter)->val << '\n';
                    }
                }

                std::deque<std::array<char, 3>> align;
                for (auto it1 = path.begin(), it2 = std::next(it1); it2 != std::prev(path.end()); ++it1, ++it2)
                {
                    if ((*it2)->w > (*it1)->w && ((*it2)->s > (*it1)->s || ((*it2)->n == (*it1)->n && (*it2)->n >= 0 && globals[(*it2)->n] && (*it2)->lambda > (*it1)->lambda)))
                    {
                        char c;
                        if (!globals[(*it2)->n])
                            c = locals[(*it2)->n]->nameseq.seq[(*it2)->s - 1];
                        else
                        {
                            auto &browheel = globals[(*it2)->n]->browheel;
                            c = BroWheel::int2base[browheel.sequence(browheel.sequence.size() - 1 - (*it2)->s)];
                        }
                        align.emplace_back(std::array<char, 3>{O[(*it2)->w - 1], '|', c});
                    }
                    else
                    {
                        for (int w = (*it1)->w; w < (*it2)->w; ++w)
                            align.emplace_back(std::array<char, 3>{O[w], ' ', '-'});
                        int n = std::max((*it1)->n, (*it2)->n);
                        if (n >= 0 && std::min((*it1)->n, (*it2)->n) < 0 && !locals[n])
                            continue;
                        if (n >= 0)
                        {
                            if (locals[n])
                            {
                                int sup = (*it2)->n < 0 ? locals[n]->nameseq.seq.size() : (*it2)->s;
                                for (int s = (*it1)->s; s < sup; ++s)
                                    align.emplace_back(std::array<char, 3>{'-', ' ', locals[n]->nameseq.seq[s]});
                            }
                            else
                            {
                                auto &browheel = globals[n]->browheel;
                                for (size_t s = (*it2)->s + (*it1)->lambda - (*it2)->lambda + 1; s <= (*it2)->s; ++s)
                                    align.emplace_back(std::array<char, 3>{'-', ' ', BroWheel::int2base[browheel.sequence(browheel.sequence.size() - 1 - s)]});
                            }
                        }
                    }
                }
                fout_align << Oname << '\t' << track_num << '\n';
                for (int i = 0; i < 3; ++i)
                {
                    for (auto &ali : align)
                        fout_align << ali[i];
                    fout_align << '\n';
                }
                Dot *dot_parent = *std::next(path.begin());
                while (path.front() == dot_parent->sources[dot_parent->s_sz - 1])
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
                        if (dot_parent->sources[i] == path.front())
                        {
                            path.pop_front();
                            path.push_front(path.front()->sources[i + 1]);
                            break;
                        }
            }
            else
                path.push_front(path.front()->sources[0]);
        }
    }

    void extract(std::string run_name, size_t max_extract)
    {
        std::ifstream fin(run_name + ".trk");
        std::ofstream fout(run_name + ".ext");
        std::string Oname;
        size_t track_num;
        size_t track_length;
        while (true)
        {
            fin >> Oname >> track_num >> track_length;
            while (fin.good() && track_num > max_extract)
            {
                for (size_t i = 0; i < track_length; ++i)
                    std::getline(fin, Oname);
                fin >> Oname >> track_num >> track_length;
            }
            if (!fin.good())
                break;
            fout << Oname << '\t' << track_num << '\n';
            if (track_length == 0)
                continue;
            else
            {
                std::string name;
                int w;
                size_t s;
                double val;
                fin >> name >> s >> w >> val;
                for (size_t i = 1; i < track_length; ++i)
                {
                    std::string namep(std::move(name));
                    int wp = w;
                    size_t sp = s;
                    double valp = val;
                    fin >> name >> s >> w >> val;
                    if (namep != name)
                    {
                        if (node_names.find(namep) != node_names.end())
                            fout << name << '\t' << s << '\t' << w << '\t' << val << '\n';
                        else
                            fout << namep << '\t' << sp << '\t' << wp << '\t' << valp << '\n';
                    }
                }
            }
        }
    }
};

#endif