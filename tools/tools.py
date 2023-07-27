import bioframe, random, Bio.Seq, pandas, subprocess, numpy, pysam, re, os, gzip, collections, more_itertools, string, re, Bio.Align

# BOWTIE2
# BOWTIE2 report up to -k default:1 concordant paired alignments with tag YT:Z:CP
# If no corcondant pair is found, it align each read, respectively, and report up to -k alignment.
# If both reads are uniquely mapped, the read pair is called discordant with tag YT:Z:DP.
# If one read is multiply-mapped, each map of both reads is reported with tag YT:Z:UP
# If say R1 is unmapped, then BOWTIE2 report R1 in a single record with tag YT:Z:UP as long as --no-unal is not set. If R2 is mapped, all of ref1/pos1, ref1N/pos1N, ref2N/pos2N are set the same as ref2/pos2. 
# Manual search of corcondant pair based on the reported YT:Z:UP records is not necessary because this has been done in the BOWITE2-internal stage of searching corcondant pair. The exception is that one manually search corcondant pair in a wider range than that specifed by -X default:500 (however, why not just set -X large enough at the first place).
# BOWTIE2 documents that specifying -k/-a will change the internal parameters swMateImmediately (whether to look for mate immediately) reportImmediately (whether to report hits immediately to msink) from false to true. Nevertheless, the fact is that swMateImmediately/reportImmediately are hard-coded as true/true.

# HISAT2
# HISAT2 defines GenomeHit<index_t>::compatibleWith(...) with the following requirements: they are on the same strand and chromosome; the hit comes first in read also comes first in reference; no hit contains the other.
# For hits satisfy GenomeHit<index_t>::compatibleWith(...) but overlap in read or reference, HISAT2 merges the hits into non-overlapping linear ones, which can be reported by a single sam record.
# Even in spliced mode (not use --no-spliced-alignment), HISAT2 still defines concordant/discordant or UP alignments recorded in YT tag for paired reads. However, -X and -I only work with -no-spliced-alignment.
# HISAT2 has a typical option --secondary. If specified, HISAT2 starts to report maps with suboptimal scores.

# STAR
# STAR reports both NH:i:<int> and HI:i:<int> tag, which can be used to group supplementary records
# STAR does not report the order of the supplementary records
# STAR also use ch tag to record all chimeric segments with --chimOutType WithinBAM and --outSAMattributes All
# STAR restricts the upper bound of read mapping loci by --outFilterMultimapNmax default:10
# STAR reports suboptimal maps controled by --outFilterMultimapScoreRange default:1 and --chimMultimapScoreRange default:1
# STAR do not report chimeric reads by default because --chimSegmentMin default:0
# R2 in paired reads is not reported as supplementary

def is_next(record1, record2):
    ### check record2 is the next record of record1 in a chimeric alignment
    fields1, fields2 = record1.split("\t"), record2.split("\t")
    qname1, flag1, chrom1, chrom1N, pos1N = fields1[0], int(fields1[1]), fields1[2], fields1[6], fields1[7]
    qname2, flag2, chrom2, pos2 = fields2[0], int(fields2[1]), fields2[2], fields2[3]

    ### check qnames
    if qname1 != qname2:
        return False
    
    ### check flags
    # check both records belong to a multiple-records chimeric group
    if not flag1 % 2 or not flag2 % 2:
        return False
    # check both records mapped
    if (flag1 / 4) % 2 or (flag2 / 4) % 2:
        return False
    # the next of record1 should map
    if (flag1 / 8) % 2:
        return False
    # record1 is the last record if and only if record2 is the first record
    if (flag1 / 128) % 2 != (flag2 / 64) % 2:
        return False

    ### check YT tag if exists
    YTpat = re.compile("\sYT:Z:([CDPU]{2})\s")
    YT1, YT2 = YTpat.search(record1), YTpat.search(record2)
    if YT1 and YT1.group(1) != "CP" or YT2 and YT2.group(1) != "CP":
        return False

    ### check chromsome 
    if chrom1N == "=":
        if chrom1 != chrom2:
            return False
    else:
        if chrom1N != chrom2:
            return False
    
    ### check position
    if pos1N != pos2:
        return False
    
    ### check strand
    if (flag1 / 32) % 2 != (flag2 / 16) % 2:
        return False
    
    ### check HI tag
    HIpat = re.compile("\sHI:i:(\d+)\s")
    HI1, HI2 = HIpat.search(record1), HIpat.search(record2)
    if HI1 and HI2 and HI1.group(1) != HI2.group(1):
        return False
    
    ### check AS and YS tags
    ASpat = re.compile("\sAS:i:(-?\d+)\s")
    YSpat = re.compile("\sYS:i:(-?\d+)\s")
    YS1, AS2 = YSpat.search(record1), ASpat.search(record2)
    if YS1 and AS2 and YS1.group(1) != AS2.group(1):
        return False
    
    return True
                        
def group_records(records, software):
    ### group records into chimeric alignments denoted by HI tag
    software = software.lower()
    if software in ["star", "hisat2"]:
        NH = int(re.search("\sNH:i:(\d+)\s", records[0]).group(1))
        HI = 0
        HIpat = re.compile("\sHI:i:(\d+)\s")
    i = 0
    while i < len(records):
        firstfound = False
        while not firstfound:
            if software in ["star", "hisat2"]:
                HI += 1
            for j in range(i, len(records)):
                flag = int(records[j].split('\t', 2)[1])
                if (flag / 64) % 2:
                    if software in ["star", "hisat2"]:
                        HIR = int(re.search("\sNH:i:(\d+)\s", records[j]).group(1))
                    if software not in ["star", "hisat2"] or HIR == HI:
                        firstfound = True
                        records[i], records[j] = records[j], records[i]
                        break
            if software not in ["star", "hisat2"] or HI == NH:
                break
        if not firstfound:
            raise Exception("the first record cannot be found")
        j = i
        flag = int(records[j].split('\t', 2)[1])
        while j + 1 < len(records) and not (flag / 128) % 2:
            for k in range(j + 1, len(records)):
                if is_next(records[j], records[k]):
                    j += 1
                    records[j], records[k] = records[k], records[j]
                    break
            flag = int(records[j].split('\t', 2)[1])
        # alignment with single record generally does not set the first and last record, so j > i is necessary
        if (j > i and not is_next(records[j], records[i])):
            raise Exception("the next of the last record is not the first record")
        if software not in ["star", "hisat2"]:
            for k in range(i, j + 1):
                records[k] = records[k][:-1] + f"\tHI:i:{HI}\n"          
        i = j + 1
    
    return records

def R1R2_unique(records, genome, aligner, extend):
    # check each read in pair is uniquely mapped, if one is uniquely mapped and the other is non-uniquely mapped, remap the other
    breakpoint()
    chroms1, poses1, strands1, chroms2, poses2, strands2 = set(), set(), set(), set(), set(), set()
    record1, record2 = "", ""
    for record in records:
        _, flag, chrom, pos, _ = record.split('\t', 4)
        flag = int(flag)
        if (flag / 16) % 2:
            strand = '-'
        else:
            strand = "+"
        if (flag / 64) % 2:
            chroms1.add(chrom)
            poses1.add(pos)
            strands1.add(strand)
            if not record1:
                record1 = record
        if (flag / 128) % 2:
            chroms2.add(chrom)
            poses2.add(pos)
            strands2.add(strand)
            if not record2:
                record2 = record
    R1unique = False
    if not chroms1:
        R1mapped = False
    else:
        R1mapped = True
        if len(chroms1) == 1 and len(poses1) == 1 and len(strands1) == 1:
            R1unique = True
    R2unique = False
    if not chroms2:
        R2mapped = False
    else:
        R2mapped = True
        if len(chroms2) == 1 and len(poses2) == 1 and len(strands2) == 1:
            R2unique = True

    remap = False
    if R1unique and R2mapped and not R2unique:
        remap = True
        recordU, recordR = record1, record2
    elif R2unique and R1mapped and not R1unique:
        remap = True
        recordU, recordR = record2, record1
    if remap:
        recordRsp = recordR.split('\t', 10)
        flagR, seq, qual = recordRsp[1], recordRsp[9], recordRsp[10]
        if (flagR / 16) % 2:
            seq = Bio.Seq.Seq(seq).reverse_complement().__str__()
            qual = qual[::-1]
        qname, flagU, chrom, posU, _, CIGARU = recordU.split('\t', 6)
        score_min = -0.6 - 0.6 * len(seq)
        if not (flagU / 16) % 2:    
            refstart = posU - 1
        else:
            tread = pysam.libcalignedsegment.AlignedSegment()
            tread.cigarstring = CIGARU
            refstart = posU - 1 + tread.reference_length - extend
        ref = genome[chrom].ff.fetch(chrom)[refstart:refstart+extend]
        if not (flagU / 16) % 2:
            seq = Bio.Seq.Seq(seq).reverse_complement().__str__()
            qual = qual[::-1]
        alignments = aligner.align(ref, seq)
        if alignments.score >= score_min:
            qname, flagU, chrom, posU, _, CIGARU, _, _, _, leftU = recordU.split('\t', 9)
            if not (flagU / 2) % 2:
                flagU += 2
            if (flagU / 8) % 2:
                flagU -= 8
            if (flagU / 16) % 2 and (flagU / 32) % 2:
                flagU -= 32
            if not (flagU / 16) % 2 and not (flagU / 16) % 2:
                flagU += 32
            if (flagU / 256) % 2:
                flagU -= 256
            leftU = re.sub("(?<=\s)YT:Z:[CDPU]{2}(?=\s)", "YT:Z:CP", leftU)
            leftU = re.sub("(?<=\s)YS:i:-?\d+(?=\s)", f"YS:i:{int(alignments.score)}", leftU)

            flagR = 1 + 2
            if (flagU / 16) % 2:
                flagR += 32
            else:
                flagR += 16
            if (flagU / 64) % 2:
                flagR += 128
            else:
                flagR += 64
            for i in range(len(alignments[0].aligned[0])):
                Trange = alignments[0].aligned[0][i]
                Qrange = alignments[0].aligned[1][i]
                if i == 0 and Qrange[0] > 0:
                    posR = refstart + Trange[0] + 1
                    CIGARR = f"{Qrange[0]}I"
                CIGARR += f"{Qrange[1] - Qrange[0]}M"
                if i < len(alignments[0].aligned[0]) - 1:
                    TrangeN = alignments[0].aligned[0][i + 1]
                    QrangeN = alignments[0].aligned[1][i + 1]
                    if TrangeN[0] > Trange[1]:
                        CIGARR += f"{TrangeN[0] - Trange[1]}D"
                    if QrangeN[0] > Qrange[1]:
                        CIGARR += f"{QrangeN[0] - Qrange[1]}I"
                else:
                    if len(seq) > Qrange[1]:
                        CIGARR += f"{len(seq) - Qrange[1]}I"
            tread.cigarstring = CIGARU
            endU = posU + tread.reference_length
            tread.cigarstring = CIGARR
            endR = posR + tread.reference_length
            TLEN = max(endU, endR) - min(posU, posR)
            if (flagU / 64) % 2:
                TLENU, TLENR = TLEN, -TLEN
            else:
                TLENR, TLENU = TLEN, -TLEN
            ASU = int(re.search("\sAS:i:(-?\d+)\s", leftU).group(1))
            recordU = f"{qname}\t{flagU}\t{chrom}\t{posU}\t255\t{CIGARU}\t{chrom}\t{posR}\t{TLENU}\t{leftU}"
            recordR = f"{qname}\t{flagR}\t{chrom}\t{posR}\t255\t{CIGARR}\t{chrom}\t{posU}\t{TLENR}\t{seq}\t{qual}\tAS:i:{int(alignments.score)}\tYS:i:{ASU}\tYT:Z:CP\n"
            if (flagU / 64) % 2:
                record1, record2 = recordU, recordR
            else:
                record1, record2 = recordR, recordU
            record1 = record1[:-1] + "\tHI:i:1\tRI:i:1\n"
            record2 = record2[:-1] + "\tHI:i:1\tRI:i:2\n"
            return record1, record2

    return None, None

def concordant_correction(samfile, outfile, oldfile, newfile):
    genome = bioframe.load_fasta("/home/ljw/hg19_with_bowtie2_index/hg19.fa")
    aligner = Bio.Align.PairwiseAligner()
    aligner.target_internal_open_gap_score = -5
    aligner.target_internal_extend_gap_score = -3
    aligner.target_left_open_gap_score = 0
    aligner.target_left_extend_gap_score = 0
    aligner.target_right_open_gap_score = 0
    aligner.target_right_extend_gap_score = 0
    aligner.query_internal_open_gap_score = -5
    aligner.query_internal_extend_gap_score = -3
    aligner.query_left_open_gap_score = -5
    aligner.query_left_extend_gap_score = -3
    aligner.query_right_open_gap_score = -5
    aligner.query_right_extend_gap_score = -3
    aligner.match_score = 0
    aligner.mismatch_score = -6
    with open(samfile, "r") as sf, open(outfile, "w") as of, open(oldfile, "w") as old, open(newfile, "w") as new:
        for record in sf:
            of.write(record)
            if (record.startswith("@HD")):
                break
        if record.find("SO:queryname") == -1:
            raise Exception("samfile must be sorted by queryname")
        for record in sf:
            if not record.startswith("@"):
                break
            else:
                of.write(record)
        qnameG = record.split('\t', 1)[0]
        records = [record]
        for record in sf:
            qname = record.split('\t', 1)[0]
            if qname == qnameG:
                records.append(record)
            else:
                record1, record2 = R1R2_unique(records, genome, aligner, 500)
                if not record1:
                    # for reco in group_records(records):
                    #     of.write(reco)
                    pass
                else:
                    of.write(record1)
                    of.write(record2)
                    for reco in group_records(records):
                        old.write(reco)
                    new.write(record1)
                    new.write(record2)
                qnameG = qname
                records = [record]
        record1, record2 = R1R2_unique(records, genome, aligner, 500)
        if not record1:
            # for reco in group_records(records):
            #     of.write(reco)
            pass
        else:
            of.write(record1)
            of.write(record2)
            for reco in group_records(records):
                old.write(reco)
            new.write(record1)
            new.write(record2)

concordant_correction("tools/tests/Z0-rep1_L3_710D01.bowtie2.sam", "tools/tests/Z0-rep1_L3_710D01.bowtie2.concordant.sam", "tools/tests/Z0-rep1_L3_710D01.bowtie2.old.sam", "tools/tests/Z0-rep1_L3_710D01.bowtie2.new.sam")

def sortedbed2bedgraph(bedfile, bedgraphfile, binwidth):
    chromsizes = {'chr1' : 249250621, 'chr2' : 243199373, 'chr3' : 198022430, 'chr4' : 191154276, 'chr5' : 180915260, 'chr6' : 171115067, 'chr7' : 159138663, 'chr8' : 146364022, 'chr9' : 141213431, 'chr10' : 135534747, 'chr11' : 135006516, 'chr12' : 133851895, 'chr13' : 115169878, 'chr14' : 107349540, 'chr15' : 102531392, 'chr16' : 90354753, 'chr17' : 81195210, 'chr18' : 78077248, 'chr19' : 59128983, 'chr20' : 63025520, 'chr21' : 48129895, 'chr22' : 51304566, 'chrM' : 16571, 'chrX' : 155270560, 'chrY' : 59373566}
    # get the coverage of a bedfile according to its scores, bedfile is assumed sorted
    bed = bioframe.read_table(bedfile, schema="bed")
    with open(bedgraphfile, "w") as bgf:
        bgf.write('track type=bedGraph name="zhangmo"\n')
        for chrom in chromsizes.keys():
            sbin, scores = 0, []
            bedchrom = bed[bed["chrom"] == chrom]
            for start, end, score in zip(bedchrom["start"], bedchrom["end"], bedchrom["score"]):
                sbinnew = start // binwidth
                ebinnew = end // binwidth
                for ib in range(sbinnew, ebinnew + 1):
                    scores.extend([0] * (ib - sbin - len(scores) + 1))
                    scores[ib - sbin] += score
                for ib in range(sbin + 1, sbinnew):
                    if scores[ib - sbin] != scores[0]:
                        if scores[0] > 0:
                            bgf.write(f"{chrom}\t{sbin * binwidth}\t{ib * binwidth}\t{scores[0]}\n")
                        scores = scores[ib - sbin:]
                        sbin = ib
            if not scores:
                continue
            i = 0
            for j in range(1, len(scores)):
                if scores[i] != scores[j]:
                    if scores[i] > 0:
                        bgf.write(f"{chrom}\t{(sbin + i) * binwidth}\t{(sbin + j) * binwidth}\t{scores[i]}\n")
                    i = j
            if scores[i] > 0:
                outend = min((j + 1) * binwidth, chromsizes[chrom])
                bgf.write(f"{chrom}\t{(sbin + i) * binwidth}\t{outend}\t{scores[i]}\n")


def count_fastq(fastq_file):
    with gzip.open(fastq_file, 'r') as f:
        duplines = collections.Counter(f.readlines()[1::4])
        with open(f'{fastq_file}.txt', 'w') as fw, open(f'{fastq_file}.fa', 'w') as fa:
            for i, seq_c in enumerate(duplines.most_common(), start=1):
                dec_seq = seq_c[0].decode()
                _ = fw.write(f"{dec_seq[:-1]}\t{seq_c[1]}\n")
                _ = fa.write(f'>seq{i}\n{dec_seq}')
    
    return f'{fastq_file}.txt', f'{fastq_file}.fa'

def to_sam(fasta_file):
    # save the genome size to ref2len
    ref2len = {}
    for row in bioframe.read_table('/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes', schema=['chrom', 'length']).itertuples():
        ref2len[row.chrom] = row.length
    # scan all .ext and .alg files
    files = [os.path.join(os.path.dirname(fasta_file),file) for file in os.listdir(os.path.dirname(fasta_file)) if file.endswith('.ext') and file.startswith(os.path.basename(fasta_file))]
    with open(f"{fasta_file}.sam", "w") as f:
        f.write("@HD\tVN:1.0\tSO:queryname\n")
        for file in files:
            # read all ext results
            with open(file, 'r') as ext:
                exts = ext.readlines()
            # scan alg line by line
            with open(f'{file[:-3]}alg', 'r') as alg:
                # use ei to index records (exts) in .ext file
                ei = -1
                for _, query, mid, ref in more_itertools.batched(alg,4):
                    query_segs = query.strip().split()
                    ref_segs = ref.strip().split()
                    ei += 1
                    ext_segs = exts[ei].strip().split()
                    QNAME = ext_segs[0][1:]
                    MAPIDX = int(ext_segs[1])
                    RNAMEs, strands, ref_starts, ref_ends, query_starts, query_ends, scores = [], [], [], [], [], [], []
                    for i in range(len(query_segs)):
                        ei += 1
                        ext_segs = exts[ei].strip().split()
                        ar = ext_segs[0].split(":")
                        if (ar[-1].endswith("_RC")):
                            RNAMEs.append(ar[-1][1:-3])
                            strands.append("-")
                            ref_ends.append(ref2len[RNAMEs[-1]]-int(ext_segs[1])+1)
                        else:
                            RNAMEs.append(ar[-1][1:])
                            strands.append("+")
                            ref_starts.append(int(ext_segs[1])+1)
                        query_starts.append(int(ext_segs[-2]))
                        scores.append(float(ext_segs[-1]))
                        ei += 1
                        ext_segs = exts[ei].strip().split()
                        if (ar[-1].endswith("_RC")):
                            ref_starts.append(ref2len[RNAMEs[-1]]-int(ext_segs[1])+1)
                        else:
                            ref_ends.append(int(ext_segs[1])+1)
                        query_ends.append(int(ext_segs[-2]))
                        scores[-1] = float(ext_segs[-1]) - scores[-1]
                    # set CIGARs
                    start = 0
                    CIGARs = []
                    for i in range(len(query_segs)):
                        ql = len(query_segs[i])
                        mid_seg = mid[start:start+len(query_segs[i])]
                        query_segs[i] = query_segs[i].lstrip("-")
                        mid_seg = mid_seg[len(mid_seg)-len(query_segs[i]):]
                        ref_segs[i] = ref_segs[i][len(ref_segs[i])-len(query_segs[i]):]
                        query_segs[i] = query_segs[i].rstrip("-")
                        if len(query_segs[i])<len(mid_seg):
                            mid_seg = mid_seg[:len(query_segs[i])-len(mid_seg)]
                            ref_segs[i] = ref_segs[i][:len(query_segs[i])-len(ref_segs[i])]
                        typepre="?"
                        cumsum=0
                        CIGARs.append("")
                        for j in range(len(query_segs[i])):
                            if mid_seg[j]=="|":
                                type="="
                            else:
                                if query_segs[i][j]=="-":
                                    type="D"
                                else:
                                    if ref_segs[i][j]=="-":
                                        type="I"
                                    else:
                                        type="X"
                            if type != typepre:
                                if cumsum > 0:
                                    CIGARs[-1] += f"{cumsum}{typepre}" 
                                    cumsum = 0
                                typepre=type
                            cumsum += 1
                        CIGARs[-1] += f"{cumsum}{type}"
                        start += ql+1
                    # set FLAGs
                    FLAGs = []
                    for i in range(len(query_segs)):
                        FLAGs.append(2)
                        if len(query_segs)>1:
                            FLAGs[-1] += 1
                        if strands[i]=="-":
                            FLAGs[-1] += 16
                        if strands[(i+1)%len(query_segs)]=="-":
                            FLAGs[-1] += 32
                        if i!=0:
                            FLAGs[-1] += 2048
                    # set TAGs
                    TAGs = []
                    for i in range(len(query_segs)):
                        TAGs.append(f"HI:i:{MAPIDX + 1}\tRI:i:{i + 1}")
                    # set TLENs
                    TLENs = []
                    if len(query_segs)==1:
                        TLENs.append(0)
                    else:
                        unify_chrom = True
                        for i in range(1,len(query_segs)):
                            if RNAMEs[i]!=RNAMEs[0]:
                                for j in range(len(query_segs)):
                                    TLENs.append(0)
                                unify_chrom = False
                        if unify_chrom:
                            leftmost = ref_starts[0]
                            rightmost = ref_ends[0]
                            rightidx = 0
                            for i in range(1,len(query_segs)):
                                if leftmost > ref_starts[i]:
                                    leftmost = ref_starts[i]
                                if rightmost < ref_ends[i]:
                                    rightmost = ref_ends[i]
                                    rightidx = i
                            TLEN = rightmost - leftmost
                            for i in range(len(query_segs)):
                                TLENs.append(TLEN)
                            TLENs[rightidx] = -TLEN
                    # set SEQs
                    SEQs = []
                    for i in range(len(query_segs)):
                        SEQs.append(query_segs[i].replace("-","").upper())
                        if strands[i]=="-":
                            SEQs[-1] = Bio.Seq.Seq(SEQs[-1]).reverse_complement().__str__()
                    # output
                    for i in range(len(query_segs)):
                        f.write(f'{QNAME}\t{FLAGs[i]}\t{RNAMEs[i]}\t{ref_starts[i]}\t255\t{CIGARs[i]}\t{RNAMEs[(i+1)%len(query_segs)]}\t{ref_starts[(i+1)%len(query_segs)]}\t{TLENs[i]}\t{SEQs[i]}\t*\t{TAGs[i]}\n')






def get_skip_blocks(record):
    _, _, _, start, _, cigarstring = record.split('\t', 6)
    cigars = cigarstring.split('N')
    skips = []
    for i in range(len(cigars)-1):
        skips.append(int(re.findall(r'\d+$', cigars[i])[0]))
        cigars[i] = cigars[i].rstrip(string.digits)
    tread = pysam.libcalignedsegment.AlignedSegment()
    blocks = [[int(start), 0]]
    for i in range(len(cigars)-1):
        tread.cigarstring = cigars[i]
        blocks[i][1] = blocks[i][0] + tread.reference_length
        blocks.append([blocks[i][1]+skips[i], 0])
    tread.cigarstring = cigars[-1]
    blocks[-1][1] = blocks[-1][0] + tread.reference_length
    return blocks

def SAM2BED(samfile, bedfile):
    HIpat = re.compile("\sHI:i:(\d+)\s")
    RIpat = re.compile("\sRI:i:(\d+)\s")
    ASpat = re.compile("\sAS:i:(-?\d+)\s")
    with open(samfile, "r") as sf, open(bedfile, "w") as bf:
        for record in sf:
            if (record.startswith("@HD")):
                break
        if record.find("SO:queryname") == -1:
            raise Exception("samfile must be sorted by queryname")
        for record in sf:
            if not record.startswith("@"):
                break
        bf.write("chrom\tstart\tend\tqname\tHI\tRI\tBI\tAS\tstrand\n")
        while record:
            # qname is query name, which may have multiple alignments
            # sam records with the same 1-based HI are chimeric parts of the same alignment
            # each record has unique 1-based RI
            # each record may have several linear map blocks (small indels are not relevant), block orders are indicated by 1-based BI
            # AS is the alignment score of the record. If records contains several blocks, then the first block inherits the AS, other blocks inherit 0
            qname, flag, chrom, _ = record.split('\t', 3)
            if (int(flag) / 16) % 2:
                strand = '-'
            else:
                strand = "+"
            HI = int(HIpat.search(record).group(1))
            RI = int(RIpat.search(record).group(1))
            AS = int(ASpat.search(record).group(1))
            blocks = get_skip_blocks(record)
            for BI, block in enumerate(blocks, start=1):
                if BI == 1:
                    ASB = AS
                else:
                    ASB = 0
                bf.write(f"{chrom}\t{block[0]}\t{block[1]}\t{qname}\t{HI}\t{RI}\t{BI}\t{ASB}\t{strand}\n")
            record = sf.readline()


def min_dis(regions1, regions2):
    # use dynamic programming to calculate the distant between two alignments
    score_mat = numpy.zeros([len(regions1)+1, len(regions2)+1])
    for i in range(1, len(regions1)+1):
        score_mat[i,0] = score_mat[i-1,0] + regions1[i]["end"] - regions1[i]["start"]
    for j in range(1, len(regions2)+1):
        score_mat[0,j] = score_mat[0,j-1] + regions2[j]["end"] - regions2[j]["start"]
    for i in range(1, len(regions1)+1):
        for j in range(1, len(regions2)+1):
            if regions1[i]["chrom"]!=regions2[j]["chrom"]:
                inter = 0
            else:
                startmax = max(regions1[i]["start"], regions2[j]["start"])
                endmin = min(regions1[i]["end"], regions2[j]["end"])
                inter = max(endmin-startmax, 0)
            len1 = regions1[i]["end"] - regions1[i]["start"]
            len2 = regions2[j]["end"] - regions2[j]["start"]
            score_mat[i,j] = max(score_mat[i-1,j]+len1, score_mat[i,j-1]+len2, score_mat[i-1,j-1]+len1+len2-inter)
    return score_mat[-1,-1]

def generate_random_seq(outfile="random.seq", mutfasta="mut.fasta", seqnum=1000, seqlen=150, overlen=20):
    # generate seqnum sequences of length seqlen with random overhangles of length overlen appended at both sides
    chrom_set = [f"chr{i}" for i in range(1,22)] + ["chrX", "chrY"]
    genome = bioframe.load_fasta("/home/ljw/hg19_with_bowtie2_index/hg19.fa")
    chrom_sizes = bioframe.read_table("/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes", schema=["chrom", "size"])
    chrom_sizes = chrom_sizes[[chrom in chrom_set for chrom in chrom_sizes["chrom"]]].reset_index(drop=True)
    total = sum(chrom_sizes["size"])
    chromseqs = {}
    for chrom in chrom_sizes["chrom"]:
        chromseqs[chrom] = genome[chrom].ff.fetch(chrom)
    chroms, starts, ends, names, strands, seqs = [], [], [], [], [], []
    for i in range(seqnum):
        rand = random.uniform(0, total)
        for chrom, size in zip(chrom_sizes["chrom"], chrom_sizes["size"]):
            rand -= size   
            if rand < 0:
                break
        chroms.append(chrom)
        while True:
            pos = random.randint(0,len(chromseqs[chrom])-seqlen)
            seq = chromseqs[chrom][pos:pos+seqlen]
            if seq != seq.upper() or "N" in seq:
                continue
            strand = "+"
            if random.uniform(0,1)<0.5:
                strand = "-"
                seq = Bio.Seq.Seq(seq).reverse_complement().__str__()
            break
        starts.append(pos)
        ends.append(pos+seqlen)
        names.append(f">seq{i}")
        strands.append(strand)
        seqs.append(seq)

    df_seqs = pandas.DataFrame({"chrom" : chroms, "start" : starts, "end" : ends, "name" : names, "strand" : strands, "seq" : seqs})

    muts = []
    for name, seq in zip(df_seqs["name"], df_seqs["seq"]):
        with open("tmp.fasta", "w") as fw:
            fw.write(f"{name}\n{seq}\n")
        subprocess.check_output(f"snpmutator -n 1 -s 10 -i 10 -d 10 -p {len(seq)} -F ./ tmp.fasta", shell=True)
        key, mut = list(bioframe.load_fasta("tmp_mutated_1.fasta").items())[0]
        overhangle = [random.choice(['A','T','C','G']) for i in range(2 * overlen)]
        muts.append("".join(overhangle[:overlen]) + mut.ff.fetch(key).lstrip(")") + "".join(overhangle[overlen:]))

    df_seqs["mut"] = muts

    df_seqs.to_csv(outfile, sep='\t', header=True, index=False)

    with open(mutfasta, "w") as fw:
        for name, mut in zip(df_seqs["name"], df_seqs["mut"]):
            fw.write(f"{name}\n{mut}\n")