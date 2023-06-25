import Bio.Seq, bioframe, numpy, scipy

poly = scipy.interpolate.lagrange(numpy.array([78,65,67,71,84]), numpy.arange(2,7))
numpy.round(poly(numpy.array([78,65,67,71,84]))).astype(numpy.uint8)
# fasta_file = "/home/ljw/hg19_with_bowtie2_index/hg19.fa"
fasta_file = "/home/ljw/hg19_with_bowtie2_index/test.fa"
genome = bioframe.load_fasta(fasta_file)
cumlen = 0
with open(f"{fasta_file}.for.rev.inv.concat", "wb") as fw, open(f"{fasta_file}.cumlen", "w") as fw2:
    for chr, value in list(genome.items())[::-1]:
        fw2.write(f">{chr}_RC\t{cumlen}\n")
        seq = value.ff.fetch(chr)[::-1].upper()
        cumlen += len(seq) + 1
        fw2.write(f">{chr}\t{cumlen}\n")
        cumlen += len(seq) + 1
        fw.write(numpy.round(poly(numpy.frombuffer(Bio.Seq.Seq(seq).reverse_complement().__str__().encode(), dtype=numpy.uint8))).astype(numpy.uint8).tobytes())
        fw.write(numpy.uint8(1).tobytes())
        fw.write(numpy.round(poly(numpy.frombuffer(seq.encode(), dtype=numpy.uint8))).astype(numpy.uint8).tobytes())
        fw.write(numpy.uint8(1).tobytes())
    fw.write(numpy.uint8(0).tobytes())