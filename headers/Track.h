#ifndef TRACK_H
#define TRACK_H

#include "Align.h"

struct Track : Memory, Graph
{
    std::vector<EdgeLocal *> locals;
    std::vector<EdgeGlobal *> globals;

    Track(int argc, char **argv, std::map<std::string, std::string> &file2seq, std::map<std::string, BroWheel> &file2browheel, size_t chunk_sz_) : Graph(argc, argv, file2seq, file2browheel)
    {
        Initial(chunk_sz_);
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

    void ReadTrack(std::vector<std::string> mg_files, std::string read_file, std::string file_track, std::string file_align, size_t max_seq, size_t max_track)
    {
        std::vector<std::ifstream> fins;
        for (auto &file : mg_files)
            fins.emplace_back(file, std::ios::binary);
        std::ifstream fin_seq(read_file);
        std::ofstream fout_track(file_track);
        std::ofstream fout_align(file_align);
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
        fin_seq >> std::ws;
        while (fin_seq.good() && ++seq_num <= max_seq)
        {
            std::string Oname, O;
            std::getline(fin_seq, Oname);
            std::getline(fin_seq, O);
            if ((read_file.size() > 3 && read_file.substr(read_file.size() - 3, 3) == ".fq") || (read_file.size() > 6 && read_file.substr(read_file.size() - 6, 6) == ".fastq"))
            {
                std::string tmp;
                std::getline(fin_seq, tmp);
                std::getline(fin_seq, tmp);
            }
            fin_seq >> std::ws;
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

                    std::list<Dot *> path;
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
    }

    void BackTrack(std::ofstream &fout_track, std::ofstream &fout_align, std::list<Dot *> &path, size_t max_track, std::string &Oname, std::string &O)
    {
        size_t track_num = 0;
        while (!path.empty() && track_num < max_track)
        {
            if (path.front()->s_sz == 0)
            {
                fout_track << Oname << '\t' << ++track_num << '\t' << path.size() << '\n';
                for (auto iter = path.begin(); iter != path.end();)
                {
                    if ((*iter)->n >= 0 && globals[(*iter)->n])
                    {
                        int n = (*iter)->n;
                        auto iter_s = iter;
                        do
                        {
                            ++iter;
                        } while (iter != path.end() && (*iter)->n == n);
                        int lambda = (*std::prev(iter))->lambda;
                        size_t s = (*std::prev(iter))->s;
                        do
                        {
                            fout_track << (*iter_s)->n << '\t' << s + (*iter_s)->lambda - lambda << '\t' << (*iter_s)->w << '\t' << (*iter_s)->val << '\n';
                            ++iter_s;
                        } while (iter_s != iter);
                    }
                    else
                    {
                        fout_track << (*iter)->n << '\t' << (*iter)->s << '\t' << (*iter)->w << '\t' << (*iter)->val << '\n';
                        ++iter;
                    }
                }
                std::list<std::array<char, 3>> align;
                for (auto it1 = path.begin(), it2 = std::next(it1); it2 != std::prev(path.end()); ++it1, ++it2)
                {
                    if ((*it2)->w > (*it1)->w && ((*it2)->s > (*it1)->s || ((*it2)->n == (*it1)->n && (*it2)->n >= 0 && globals[(*it2)->n] && (*it2)->lambda > (*it1)->lambda)))
                    {
                        char c;
                        if (!globals[(*it2)->n])
                            c = locals[(*it2)->n]->seq[(*it2)->s - 1];
                        else
                        {
                            auto &browheel = globals[(*it2)->n]->browheel;
                            c = BroWheel::int2base[browheel.sequence(browheel.start_rev((*it2)->s, 0))];
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
                                int sup = (*it2)->n < 0 ? locals[n]->seq.size() : (*it2)->s;
                                for (int s = (*it1)->s; s < sup; ++s)
                                    align.emplace_back(std::array<char, 3>{'-', ' ', locals[n]->seq[s]});
                            }
                            else
                            {
                                auto &browheel = globals[n]->browheel;
                                for (size_t s = (*it2)->s + (*it1)->lambda - (*it2)->lambda + 1; s <= (*it2)->s; ++s)
                                    align.emplace_back(std::array<char, 3>{'-', ' ', BroWheel::int2base[browheel.sequence(browheel.start_rev(s, 0))]});
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

    void extract(std::string file_track, std::string file_extract, size_t max_extract)
    {
        std::ifstream fin(file_track);
        std::ofstream fout(file_extract);
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
                int n, w;
                size_t s;
                double val;
                fin >> n >> s >> w >> val;
                for (size_t i = 1; i < track_length; ++i)
                {
                    int np = n;
                    int wp = w;
                    size_t sp = s;
                    fin >> n >> s >> w >> val;
                    if (np < 0 && n >= 0)
                        fout << n << '\t' << s << '\t' << w << '\n';
                    else if (np >= 0 && n < 0)
                        fout << np << '\t' << sp << '\t' << wp << '\n';
                }
            }
        }
    }
};

#endif