#ifndef HELP_H
#define HELP_H
void help()
{
    std::cout << "run = str;\n"
              << "// str is the prefix of all outputs\n"
              << "read_files = {\n"
              << "\tstr;\n"
              << "\tstr;\n"
              << "\t...\n"
              << "};\n"
              << "// the input reads\n"
              << "threads_sz = int;\n"
              << "// the number of worker threads (not including the main thread) used to align the reads\n"
              << "max_round = int;\n"
              << "// the maximal possible round (start from 0)\n"
              << "max_extract = int;\n"
              << "// the maximal number of alignments to extract\n"
              << "diff_thres = int;\n"
              << "// the minimal difference in base to distinct two extracts\n"
              << "max_range = int;\n"
              << "// the maximal number of range to use from a single suffix\n"
              << "min_seg_num = int;\n"
              << "// segment number must be larger than this\n"
              << "max_seg_num = int;\n"
              << "// segment number must be smaller than this\n"
              << "max_mega = int;\n"
              << "// maximal number of end points of segments to track\n"
              << "block_size = int;\n"
              << "// the size of read batch distributed to an alignment workder\n"
              << "nodes = {\n"
              << "\t{\n"
              << "\t\tname = str;\n"
              << "\t\t// node name\n"
              << "\t\tve = double;\n"
              << "\t\t// gap-open penalty between segments before and after this node\n"
              << "\t\tue = double;\n"
              << "\t\t// gap-extend penalty between segments before and after this node\n"
              << "\t};\n"
              << "\t...\n"
              << "};\n"
              << "// node information\n"
              << "roots = {\n"
              << "\tstr;\n"
              << "\tstr;\n"
              << "\t...\n"
              << "};\n"
              << "// names of roots of the graph\n"
              << "targets = {\n"
              << "\tstr;\n"
              << "\tstr;\n"
              << "\t...\n"
              << "};\n"
              << "// names of targets of the graph\n"
              << "locals = {\n"
              << "\t{\n"
              << "\t\tgamma = default;\n"
              << "\t\tgamma = {\n"
              << "\t\t\tKK;\n"
              << "\t\t\tKL;\n"
              << "\t\t\tKN;\n"
              << "\t\t\tKA;\n"
              << "\t\t\tKC;\n"
              << "\t\t\tKT;\n"
              << "\t\t\tKG;\n"
              << "\t\t\tLK;\n"
              << "\t\t\tLL;\n"
              << "\t\t\tLN;\n"
              << "\t\t\tLA;\n"
              << "\t\t\tLC;\n"
              << "\t\t\tLT;\n"
              << "\t\t\tLG;\n"
              << "\t\t\tNK;\n"
              << "\t\t\tNL;\n"
              << "\t\t\tNN;\n"
              << "\t\t\tNA;\n"
              << "\t\t\tNC;\n"
              << "\t\t\tNT;\n"
              << "\t\t\tNG;\n"
              << "\t\t\tAK;\n"
              << "\t\t\tAL;\n"
              << "\t\t\tAN;\n"
              << "\t\t\tAA;\n"
              << "\t\t\tAC;\n"
              << "\t\t\tAT;\n"
              << "\t\t\tAG;\n"
              << "\t\t\tCK;\n"
              << "\t\t\tCL;\n"
              << "\t\t\tCN;\n"
              << "\t\t\tCA;\n"
              << "\t\t\tCC;\n"
              << "\t\t\tCT;\n"
              << "\t\t\tCG;\n"
              << "\t\t\tTK;\n"
              << "\t\t\tTL;\n"
              << "\t\t\tTN;\n"
              << "\t\t\tTA;\n"
              << "\t\t\tTC;\n"
              << "\t\t\tTT;\n"
              << "\t\t\tTG;\n"
              << "\t\t\tGK;\n"
              << "\t\t\tGL;\n"
              << "\t\t\tGN;\n"
              << "\t\t\tGA;\n"
              << "\t\t\tGC;\n"
              << "\t\t\tGT;\n"
              << "\t\t\tGG;\n"
              << "\t\t};\n"
              << "\t\t// use the defauly blast (mis)match score parameters or specify by hand\n"
              << "\t\tname = str;\n"
              << "\t\t// path to the global reference\n"
              << "\t\tve = double;\n"
              << "\t\t// horizonal gap-open penalty, usually negative\n"
              << "\t\tue = double;\n"
              << "\t\t// horizonal gap-extend penalty, usually negative\n"
              << "\t\tvf = double;\n"
              << "\t\t// vertical gap-open penalty, usually negative\n"
              << "\t\tuf = double;\n"
              << "\t\t// vertical gap-extend penalty, usually negative\n"
              << "\t\tT = double;\n"
              << "\t\t// basic penalty, usually negative\n"
              << "\t\tdt = double;\n"
              << "\t\t// increment of basic penalty in the loop strategy, usually negative\n"
              << "\t\tmin_score = double;\n"
              << "\t\t// minimal-score attribute when this edge is in a circuit\n"
              << "\t\ttail = str;\n"
              << "\t\t// tail node name\n"
              << "\t\thead = str;\n"
              << "\t\t// head node name\n"
              << "\t\tvfp = double;\n"
              << "\t\t// gap-open penalty at the beginning of reference\n"
              << "\t\tufp = double;\n"
              << "\t\t// gap-extend penalty at the beginning of reference\n"
              << "\t\tvfm = double;\n"
              << "\t\t// gap-open penalty at the end of reference\n"
              << "\t\tufm = double;\n"
              << "\t\t// gap-extend penalty at the end of reference\n"
              << "\t};\n"
              << "\t...\n"
              << "};\n"
              << "// local-edge information\n"
              << "globals = {\n"
              << "\t{\n"
              << "\t\tgamma = default;\n"
              << "\t\tgamma = {\n"
              << "\t\t\t// similar as locals;\n"
              << "\t\t};\n"
              << "\t\t// use the defauly blast (mis)match score parameters or specify by hand\n"
              << "\t\tname = str;\n"
              << "\t\t// path to the global reference\n"
              << "\t\tve = double;\n"
              << "\t\t// horizonal gap-open penalty, usually negative\n"
              << "\t\tue = double;\n"
              << "\t\t// horizonal gap-extend penalty, usually negative\n"
              << "\t\tvf = double;\n"
              << "\t\t// vertical gap-open penalty, usually negative\n"
              << "\t\tuf = double;\n"
              << "\t\t// vertical gap-extend penalty, usually negative\n"
              << "\t\tT = double;\n"
              << "\t\t// basic penalty, usually negative\n"
              << "\t\tdt = double;\n"
              << "\t\t// increment of basic penalty in the loop strategy, usually negative\n"
              << "\t\tmin_score = double;\n"
              << "\t\t// minimal-score attribute when this edge is in a circuit\n"
              << "\t\ttail = str;\n"
              << "\t\t// tail node name\n"
              << "\t\thead = str;\n"
              << "\t\t// head node name\n"
              << "\t\treverse_complement = true/false;\n"
              << "\t\t// whether also map to the reverse-complement of the referecen\n"
              << "\t\t// this influences indexing reference, so will only take effect after reindexing\n"
              << "\t};\n"
              << "\t...\n"
              << "};\n"
              << "// global-edge information\n";
}
#endif