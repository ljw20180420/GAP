#ifndef RANDOMREADS_H
#define RANDOMREADS_H
#include <cstdlib>
#include <utility>
#include <set>
#include <iterator>
#include "Track.h"

std::tuple<std::string, std::string, std::string> random_mut(std::string &str, double indel_rate, double mut_rate, int head_in, int tail_in)
{
    std::string mut;
    int64_t row = 0;
    while (row <= str.size())
    {
        double rv = double(rand()) / RAND_MAX;
        if (rv < indel_rate)
            mut.push_back(BroWheel::int2base[rand() % 4 + 3]);
        else
        {
            ++row;
            if (row <= str.size() && rv >= 2 * indel_rate)
            {
                if (rv < 2 * indel_rate + mut_rate)
                    mut.push_back(BroWheel::int2base[rand() % 4 + 3]);
                else
                    mut.push_back(str[row - 1]);
            }
        }
    }
    std::string head_str, tail_str;
    head_in = rand() % head_in;
    tail_in = rand() % tail_in;
    for (int i = 0; i < head_in; ++i)
        head_str.push_back(BroWheel::int2base[rand() % 4 + 3]);
    for (int i = 0; i < tail_in; ++i)
        tail_str.push_back(BroWheel::int2base[rand() % 4 + 3]);
    return std::make_tuple(head_str, mut, tail_str);
}

void random_seq(Graph &graph, std::ofstream &fout_read, std::ofstream &fout_truth, int aseqlb, int asequb, double indel_rate, double mut_rate, int head_in, int tail_in, double apro)
{
    std::string seq;
    Node *node = graph.roots[rand() % graph.roots.size()];
    while (node->out_local_crosses.size() + node->out_local_circuits.size() + node->out_global_crosses.size() + node->out_global_circuits.size() == 0)
        node = graph.roots[rand() % graph.roots.size()];
    bool reach_target = false;
    do
    {
        for (Node *target : graph.targets)
            if (node->name == target->name)
            {
                if ((node->out_local_crosses.size() + node->out_local_circuits.size() + node->out_global_crosses.size() + node->out_global_circuits.size() == 0) || double(rand()) / RAND_MAX < apro)
                    reach_target = true;
                break;
            }
        if (reach_target)
            break;
        int64_t rv = rand() % (node->out_local_crosses.size() + node->out_local_circuits.size() + node->out_global_crosses.size() + node->out_global_circuits.size());
        int64_t len, start;
        std::string str, head_str, mut, tail_str;
        if (rv < node->out_local_crosses.size() + node->out_local_circuits.size())
        {
            EdgeLocal *local;
            if (rv < node->out_local_crosses.size())
                local = node->out_local_crosses[rv];
            else
                local = node->out_local_circuits[rv - node->out_local_crosses.size()];
            len = aseqlb + rand() % (asequb - aseqlb);
            start = rand() % (local->nameseq.seq.size() - len + 1);
            str = local->nameseq.seq.substr(start, len);
            std::tie(head_str, mut, tail_str) = random_mut(str, indel_rate, mut_rate, head_in, tail_in);
            seq += head_str;
            fout_truth << local->name << ':' << local->nameseq.name << '\t' << start << '\t' << seq.size() << '\n';
            seq += mut;
            fout_truth << local->name << ':' << local->nameseq.name << '\t' << start + len << '\t' << seq.size() << '\n';
            seq += tail_str;
            node = local->head;
        }
        else
        {
            EdgeGlobal *global;
            if (rv < node->out_local_crosses.size() + node->out_local_circuits.size() + node->out_global_crosses.size())
                global = node->out_global_crosses[rv - node->out_local_crosses.size() - node->out_local_circuits.size()];
            else
                global = node->out_global_circuits[rv - node->out_local_crosses.size() - node->out_local_circuits.size() - node->out_global_crosses.size()];
            len = aseqlb + rand() % (asequb - aseqlb);
            start = rand() % (global->browheel.sequence.size() - len);
            for (int64_t i = 1; i <= len; ++i)
                str.push_back(BroWheel::int2base[global->browheel.sequence(global->browheel.sequence.size() - 1 - start - i)]);
            std::tie(head_str, mut, tail_str) = random_mut(str, indel_rate, mut_rate, head_in, tail_in);
            seq += head_str;
            fout_truth << global->name << ':' << global->browheel.name_cumlen.front().first << '\t' << start << '\t' << seq.size() << '\n';
            seq += mut;
            fout_truth << global->name << ':' << global->browheel.name_cumlen.front().first << '\t' << start + len << '\t' << seq.size() << '\n';
            seq += tail_str;
            node = global->head;
        }
    } while (true);
    fout_read << seq << '\n';
}

void random_DG(int n_sz, int r_sz, int t_sz, int max_e_sz, std::string argfile, bool acyclic, double gpro, double rpro, double nve, double nue, double ve, double ue, double vf, double uf, double T, double vfp, double ufp, double vfm, double ufm, double mat, double mis, std::string local_file, std::string global_file, int lseqlb, int lsequb, int gseqlb, int gsequb, int aseqlb, int asequb, int seq_num, std::string read_file, std::string truth_file, double indel_rate, double mut_rate, int head_in, int tail_in, double apro)
{
    srand(1);
    std::vector<std::vector<std::string>> args;
    args.emplace_back(std::vector<std::string>{"--nodes"});
    for (int i = 0; i < n_sz; ++i)
        args.emplace_back(std::vector<std::string>{"node" + std::to_string(i), std::to_string(nve), std::to_string(nue)});
    args.emplace_back(std::vector<std::string>{"--roots"});
    args.emplace_back(std::vector<std::string>{"node0"});
    std::vector<int64_t> indices(n_sz - 1);
    std::iota(indices.begin(), indices.end(), 1);
    std::random_shuffle(indices.begin(), indices.end());
    for (int i = 0; i < r_sz - 1; ++i)
        args.emplace_back(std::vector<std::string>{"node" + std::to_string(indices[i])});
    std::random_shuffle(indices.begin(), indices.end());
    args.emplace_back(std::vector<std::string>{"--targets"});
    args.emplace_back(std::vector<std::string>{"node" + std::to_string(n_sz - 1)});
    for (int i = 0; i < t_sz - 1; ++i)
        args.emplace_back(std::vector<std::string>{"node" + std::to_string(indices[i])});
    std::vector<std::vector<std::string>> args_local;
    std::vector<std::vector<std::string>> args_global;
    std::map<std::string, NameSeq> file2seq;
    std::map<std::string, BroWheel> file2browheel;
    int argc;
    char **argv;
    bool at_least_two_sccs = false, average_four_edge_types;
    do
    {
        args_local.clear();
        args_global.clear();
        file2seq.clear();
        file2browheel.clear();
        while (true)
        {
            int h, t;
            do
            {
                h = rand() % n_sz;
                t = rand() % n_sz;
            } while (acyclic && h <= t);
            bool rep = false;
            for (auto &arg : args_local)
                if (arg.back() == "node" + std::to_string(h) && arg[arg.size() - 2] == "node" + std::to_string(t))
                {
                    rep = true;
                    break;
                }
            if (!rep)
                for (auto &arg : args_global)
                    if (arg.back() == "node" + std::to_string(h) && arg[arg.size() - 2] == "node" + std::to_string(t))
                    {
                        rep = true;
                        break;
                    }
            if (!rep || double(rand()) / RAND_MAX < rpro)
            {
                int en = args_local.size() + args_global.size();
                std::string mat_s = std::to_string(mat), mis_s = std::to_string(mis);
                if (double(rand()) / RAND_MAX < gpro)
                {
                    args_global.emplace_back(std::vector<std::string>{
                        "-inf", "-inf", "-inf", "-inf", "-inf", "-inf", "-inf",
                        "-inf", "-inf", "-inf", "-inf", "-inf", "-inf", "-inf",
                        "-inf", "-inf", "-inf", "-inf", "-inf", "-inf", "-inf",
                        "-inf", "-inf", "-inf", mat_s, mis_s, mis_s, mis_s,
                        "-inf", "-inf", "-inf", mis_s, mat_s, mis_s, mis_s,
                        "-inf", "-inf", "-inf", mis_s, mis_s, mat_s, mis_s,
                        "-inf", "-inf", "-inf", mis_s, mis_s, mis_s, mat_s,
                        global_file + std::to_string(en), std::to_string(ve), std::to_string(ue), std::to_string(vf), std::to_string(uf), std::to_string(T), "node" + std::to_string(t), "node" + std::to_string(h)});
                }
                else
                {
                    args_local.emplace_back(std::vector<std::string>{
                        "-inf", "-inf", "-inf", "-inf", "-inf", "-inf", "-inf",
                        "-inf", "-inf", "-inf", "-inf", "-inf", "-inf", "-inf",
                        "-inf", "-inf", "-inf", "-inf", "-inf", "-inf", "-inf",
                        "-inf", "-inf", "-inf", mat_s, mis_s, mis_s, mis_s,
                        "-inf", "-inf", "-inf", mis_s, mat_s, mis_s, mis_s,
                        "-inf", "-inf", "-inf", mis_s, mis_s, mat_s, mis_s,
                        "-inf", "-inf", "-inf", mis_s, mis_s, mis_s, mat_s,
                        local_file + std::to_string(en), std::to_string(ve), std::to_string(ue), std::to_string(vf), std::to_string(uf), std::to_string(T), std::to_string(vfp), std::to_string(ufp), std::to_string(vfm), std::to_string(ufm), "node" + std::to_string(t), "node" + std::to_string(h)});
                }
            }
            else
                continue;
            argc = 3;
            for (auto &aa : {args, args_local, args_global})
                for (auto &a : aa)
                    argc += a.size();
            argv = new char *[argc];
            int ii = 0;
            for (auto &arg : args)
                for (auto &ar : arg)
                    argv[++ii] = (char *)ar.c_str();
            argv[++ii] = (char *)"--locals";
            for (auto &arg : args_local)
                for (auto &ar : arg)
                    argv[++ii] = (char *)ar.c_str();
            argv[++ii] = (char *)"--globals";
            for (auto &arg : args_global)
                for (auto &ar : arg)
                    argv[++ii] = (char *)ar.c_str();
            Graph graph(argc, argv, file2seq, file2browheel);
            bool not_connect = false;
            for (auto &node : graph.nodes)
            {
                if (node.in_global_circuits.size() + node.in_global_crosses.size() + node.in_local_circuits.size() + node.in_local_crosses.size() + node.out_global_circuits.size() + node.out_global_crosses.size() + node.out_local_circuits.size() + node.out_local_crosses.size() == 0)
                {
                    not_connect = true;
                }
            }
            if (!not_connect)
            {
                if (graph.sccs.size() > 1)
                    at_least_two_sccs = true;
                int64_t up_single_type = (args_local.size() + args_global.size() + 2) / 3;
                average_four_edge_types = graph.local_crosses.size() > 0 && graph.local_crosses.size() <= up_single_type && graph.local_circuits.size() > 0 && graph.local_circuits.size() <= up_single_type && graph.global_crosses.size() > 0 && graph.global_crosses.size() <= up_single_type && graph.global_circuits.size() > 0 && graph.global_circuits.size() <= up_single_type;
                break;
            }
            delete[] argv;
        }
    } while (!at_least_two_sccs || args_local.size() + args_global.size() > max_e_sz || !average_four_edge_types);

    for (auto &pair : file2seq)
    {
        int64_t seqlen = lseqlb + rand() % (lsequb - lseqlb);
        std::ofstream fout(pair.first);
        fout << '>' << pair.first << '\n';
        for (int64_t j = 0; j < seqlen; ++j)
            fout << BroWheel::int2base[rand() % 4 + 3];
        fout << '\n';
        fout.close();
        pair.second.readin(pair.first);
    }
    for (auto &pair : file2browheel)
    {
        int64_t seqlen = gseqlb + rand() % (gsequb - gseqlb);
        std::ofstream fout(pair.first);
        fout << '>' << pair.first << '\n';
        for (int64_t j = 0; j < seqlen; ++j)
            fout << BroWheel::int2base[rand() % 4 + 3];
        fout << '\n';
        fout.close();
        pair.second.readin(pair.first, true, false);
    }
    Graph graph(argc, argv, file2seq, file2browheel);
    delete[] argv;
    std::ofstream fout(argfile);
    fout << "\"--read_files\", "
         << "\"" << read_file << "\",\n";
    for (auto &arg : args)
    {
        for (auto &ar : arg)
            fout << "\"" << ar << "\", ";
        fout << '\n';
    }
    fout << "\"--locals\",\n";
    for (auto &arg : args_local)
    {
        for (auto &ar : arg)
            fout << "\"" << ar << "\", ";
        fout << '\n';
    }
    fout << "\"--globals\",\n";
    for (auto &arg : args_global)
    {
        for (auto &ar : arg)
            fout << "\"" << ar << "\", ";
        fout << '\n';
    }
    fout.close();

    std::ofstream fout_read(read_file), fout_truth(truth_file);
    for (int sn = 1; sn <= seq_num; ++sn)
    {
        fout_read << ">seq" << sn << '\n';
        fout_truth << ">seq" << sn << '\n';
        random_seq(graph, fout_read, fout_truth, aseqlb, asequb, indel_rate, mut_rate, head_in, tail_in, apro);
    }
    fout_read.close();
    fout_truth.close();
}

#endif