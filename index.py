import Bio.Seq, bioframe, ctypes, struct



# genome = bioframe.load_fasta("/home/ljw/hg19_with_bowtie2_index/hg19.fa")
# batch = 50
# with open("/home/ljw/hg19_with_bowtie2_index/hg19.with.revcomp.fa", "w") as fw:
#     for chr, value in genome.items():
#         fw.write(f">{chr}\n")
#         seq = value.ff.fetch(chr)
#         seq = "\n".join([seq[i:min(i+50,len(seq))] for i in range(0,len(seq),batch)])
#         fw.write(f"{seq}\n")
#         fw.write(f">{chr}_RC\n")
#         seq = Bio.Seq.Seq(value.ff.fetch(chr)).reverse_complement().__str__()
#         seq = "\n".join([seq[i:min(i+batch,len(seq))] for i in range(0,len(seq),batch)])
#         fw.write(f"{seq}\n")
        


# with open("/home/ljw/hg19_with_bowtie2_index/hg19.fa.8.8.gsa", "rb") as f:
#     for _ in range(20):
#         struct.unpack('Q', f.read(8))

# with open("/home/ljw/hg19_with_bowtie2_index/test.fa", "r") as fr, open("/home/ljw/hg19_with_bowtie2_index/test.rev.fa", "w") as fw:
#     lines = fr.readlines()
#     lines.reverse()
#     for i in range(0,len(lines), 2):
#         seq = lines[i].rstrip('\n')
#         name = lines[i+1].rstrip('\n')
#         fw.write(f'{name}_R\n{seq[::-1]}\n')

def ca2gsu(seq):
    mydict = {"K" : "#", "L" : "$", "N" : "%", "A" : "A", "C" : "B", "T" : "C", "G" : "D"}
    nseq = ""
    for ca in seq:
        nseq += mydict[ca]
    return nseq

# def ca2ABCD(seq):
#     mydict = {"A" : "A", "C" : "B", "T" : "C", "G" : "D"}
#     nseq = ""
#     for ca in seq:
#         nseq += mydict[ca]
#     return nseq

# def num2gsu(seq):
#     mydict = {"0" : "#", "1" : "$", "2" : "N", "3" : "A", "4" : "B", "5" : "C" , "6" : "D"}
#     nseq = ""
#     for ca in seq:
#         nseq += mydict[ca]
#     return nseq

with open("/home/ljw/hg19_with_bowtie2_index/test.rev.fa", "r") as fr:
    lines = fr.readlines()
catline = ""
for i in range(1,len(lines),2):
    catline += lines[i].rstrip('\n')
    catline += 'L'
catline = catline + 'K'
catline = ca2gsu(catline)

suffices = []
for i in range(len(catline)):
    suffices.append(catline[i:] + catline[:i])
suffices.sort()

bwt = ""
for suffix in suffices:
    bwt += suffix[-1]

with open("/home/ljw/hg19_with_bowtie2_index/test.fa.my.bwt", "r") as fr:
    mybwt = fr.readline().rstrip('\n')
mybwt = ca2gsu(mybwt)

with open("/home/ljw/hg19_with_bowtie2_index/test.rev.fa", "r") as fr, open("/home/ljw/hg19_with_bowtie2_index/test.rev.gsu.fa", "w") as fw:
    lines = fr.readlines()
    for i in range(0,len(lines),2):
        fw.write(f"{lines[i]}{ca2gsu(lines[i+1][:-1])}\n")


with open("/home/ljw/hg19_with_bowtie2_index/test.fa", "r") as fr, open("/home/ljw/hg19_with_bowtie2_index/test.low.fa", "w") as fw:
    lines = fr.readlines()
    for i in range(0,len(lines),2):
        fw.write(f"{lines[i]}{lines[i+1][:-1].lower()}\n")