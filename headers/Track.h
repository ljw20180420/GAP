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
        for (auto &edge : local_crosses)
        {
            if (edge.n>=locals.size())
                locals.resize(edge.n+1,NULL);
            locals[edge.n]=&edge;
        }
        for (auto &edge : local_circuits)
        {
            if (edge.n>=locals.size())
                locals.resize(edge.n+1,NULL);
            locals[edge.n]=&edge;
        }
        for (auto &edge : global_crosses)
        {
            if (edge.n>=globals.size())
                globals.resize(edge.n+1,NULL);
            globals[edge.n]=&edge;
        }
        for (auto &edge : global_circuits)
        {
            if (edge.n>=globals.size())
                globals.resize(edge.n+1,NULL);
            globals[edge.n]=&edge;
        }
    }

    void ReadTrack(std::vector<std::string> files, std::string file_seq, std::string file_track, std::string file_align, size_t max_seq, size_t max_track)
    {
        std::vector<std::ifstream> fins;
        for (auto &file : files)
            fins.emplace_back(file, std::ios::binary);
        std::ifstream fin_seq(file_seq);
        std::ofstream fout_track(file_track);
        std::ofstream fout_align(file_align);
        size_t max_id = 0;
        std::vector<std::string> Onames;
        for (int i = 0; i < fins.size(); ++i)
        {
            size_t id;
            fins[i].seekg(-sizeof(size_t), fins[i].end);
            fins[i].read((char *)&id, sizeof(size_t));
            max_id = id > max_id ? id : max_id;
            fins[i].seekg(0, fins[i].beg);
            if (fins[i].end - fins[i].tellg() > sizeof(size_t))
            {
                int Onsz;
                fins[i].read((char *)&Onsz, sizeof(Onsz));
                Onames[i].resize(Onsz);
                fins[i].read(Onames[i].data(), sizeof(char) * Onsz);
            }
            else
                fins[i].setstate(std::ios::failbit);
        }
        std::vector<Dot> dots(max_id + 1);
        size_t seq_num=0;
        fin_seq >> std::ws;
        while (fin_seq.good() && ++seq_num<=max_seq)
        {
            std::string Oname, O;
            std::getline(fin_seq, Oname);
            std::getline(fin_seq, O);
            if ((file_seq.size() > 3 && file_seq.substr(file_seq.size() - 3, 3) == ".fq") || (file_seq.size() > 6 && file_seq.substr(file_seq.size() - 6, 6) == ".fastq"))
            {
                std::string tmp;
                std::getline(fin_seq, tmp);
                std::getline(fin_seq, tmp);
            }
            fin_seq >> std::ws;
            for (int i = 0; i < fins.size(); ++i)
                if (!fins[i].fail() && Onames[i] == Oname)
                {
                    ssize_t id = -1, uid=0;
                    do
                    {
                        ++id;
                        fins[i].read((char *)&dots[id].n, sizeof(Dot::n));
                        fins[i].read((char *)&dots[id].s, sizeof(Dot::s));
                        fins[i].read((char *)&dots[id].w, sizeof(Dot::w));
                        fins[i].read((char *)&dots[id].val, sizeof(Dot::val));
                        fins[i].read((char *)&dots[id].s_sz, sizeof(Dot::s_sz));
                        dots[id].sources = heap_alloc<Dot *>(dots[id].s_sz);
                        for (int i = 0; i < dots[id].s_sz; ++i)
                        {
                            size_t idd;
                            fins[i].read((char *)&idd, sizeof(idd));
                            dots[id].sources[i] = &dots[idd];
                            uid = idd>uid? idd : uid;
                        }
                    } while (id < uid);

                    std::list<Dot *> path;
                    path.push_front(&dots[0]);
                    size_t track_num=0;
                    BackTrack(fout_track, fout_align, path, track_num, max_track, Oname, O);

                    clear();
                    if (fins[i].end - fins[i].tellg() > sizeof(size_t))
                    {
                        int Onsz;
                        fins[i].read((char *)&Onsz, sizeof(Onsz));
                        Onames[i].resize(Onsz);
                        fins[i].read(Onames[i].data(), sizeof(char) * Onsz);
                    }
                    else
                        fins[i].setstate(std::ios::failbit);
                    break;
                }
        }
    }

    void BackTrack(std::ofstream &fout_track, std::ofstream &fout_align, std::list<Dot *> &path, size_t &track_num, size_t max_track, std::string &Oname, std::string &O)
    {
        if (path.front()->s_sz == 0)
        {
            fout_track << Oname << '\t' << ++track_num << '\t' << path.size() << '\n';
            for (auto dot : path)
                fout_track << dot->n << '\t' << dot->s << '\t' << dot->n << '\n';
            std::list<std::array<char, 3>> align;
            for (auto it1 = path.begin(), it2 = std::next(it1); it2 != path.end(); ++it1, ++it2)
            {
                if ((*it2)->s != (*it1)->s && (*it2)->w != (*it1)->w)
                {
                    char c;
                    if (globals[(*it1)->n]==NULL)
                        c=locals[(*it1)->n]->seq[(*it1)->s];
                    else
                    {
                        auto &browheel=globals[(*it1)->n]->browheel;
                        c=BroWheel::int2base[browheel.sequence(browheel.start_rev((*it1)->s,1))];
                    }
                    align.emplace_back(O[(*it1)->w],'|',c);
                }
                else
                {
                    for (int w=(*it1)->w; w<(*it2)->w; ++w)
                        align.emplace_back(O[w],' ','-');
                    int n=std::max((*it1)->n,(*it2)->n);
                    if (n>=0 && std::min((*it1)->n,(*it2)->n)<0 && locals[n]==NULL)
                        continue;
                    if (n>=0)
                    {
                        if (locals[n]!=NULL)
                        {
                            int sup= (*it2)->n<0 ? locals[n]->seq.size() : (*it2)->s;
                            for (int s=(*it1)->s; s<sup; ++s)
                                align.emplace_back('-',' ',locals[n]->seq[s]);
                        }
                        else
                        {
                            auto &browheel=globals[n]->browheel;
                            for (int s=(*it1)->s; s<(*it2)->s; ++s)
                                align.emplace_back('-',' ',BroWheel::int2base[browheel.sequence(browheel.start_rev(s,1))]);
                        }
                    }
                }
            }
            fout_align << Oname << '\t' << track_num << '\n';
            for (int i=0; i<3; ++i)
            {
                for (auto &ali : align)
                    fout_align << ali[i];
                fout_align << '\n';
            }
        }
        else
        {
            for (int i=0; i<path.front()->s_sz; ++i)
            {
                path.push_front(path.front()->sources[i]);
                BackTrack(fout_track, fout_align, path, track_num, max_track, Oname, O);
                if (track_num==max_track)
                    break;
                path.pop_front();
            }
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
            while (fin.good() && track_num>max_extract)
            {
                for (size_t i=0; i<track_length; ++i)
                    std::getline(fin, Oname);
                fin >> Oname >> track_num >> track_length;
            }
            if (!fin.good())
                break;
            fout << Oname << '\t' << track_num <<'\n';
            if (track_length==0)
                continue;
            else
            {
                int n, w;
                size_t s;
                fin >> n >> s >> w;
                for (size_t i=1; i<track_length; ++i)
                {
                    int np=n;
                    int wp=w;
                    size_t sp=s;
                    fin >> n >> s >> w;
                    if (np<0 && n>=0)
                        fout << n << '\t' << s << '\t' << w <<'\n';
                    else if (np>=0 && n<0)
                        fout << np << '\t' << sp << '\t' << wp <<'\n';
                }
            }
        }
    }
};

#endif