import Bio.Seq, bioframe, numpy, scipy

def numerize_concat_fasta(fasta_file, poly, inverse, reverse_complement, addDollar, addSharp):
    genome = bioframe.load_fasta(fasta_file)
    cumlen = 0
    with open(f"{fasta_file}.{inverse}.{reverse_complement}.{addDollar}.{addSharp}.concat", "wb") as fw, open(f"{fasta_file}.{inverse}.{reverse_complement}.{addDollar}.{addSharp}.concat.cumlen", "w") as fw2:
        if inverse:
            gitems = list(genome.items())[::-1]
        else:
            gitems = list(genome.items())
        for chr, value in gitems:
            chrs = [chr]
            if reverse_complement:
                chrs += [f"{chr}_RC"]
            if inverse:
                chrs = chrs[::-1]
            fw2.write(f">{chrs[0]}\t{cumlen}\n")
            if inverse:
                seq = value.ff.fetch(chr)[::-1].upper()
            else:
                seq = value.ff.fetch(chr).upper()
            cumlen += len(seq) + addDollar
            if reverse_complement:
                fw2.write(f">{chrs[1]}\t{cumlen}\n")
                cumlen += len(seq) + addDollar

            if inverse and reverse_complement:
                fw.write(numpy.round(poly(numpy.frombuffer(Bio.Seq.Seq(seq).reverse_complement().__str__().encode(), dtype=numpy.uint8))).astype(numpy.uint8).tobytes())
                if addDollar:
                    fw.write(numpy.uint8(1).tobytes())
            fw.write(numpy.round(poly(numpy.frombuffer(seq.encode(), dtype=numpy.uint8))).astype(numpy.uint8).tobytes())
            if addDollar:
                fw.write(numpy.uint8(1).tobytes())
            if not inverse and reverse_complement:
                fw.write(numpy.round(poly(numpy.frombuffer(Bio.Seq.Seq(seq).reverse_complement().__str__().encode(), dtype=numpy.uint8))).astype(numpy.uint8).tobytes())
                if addDollar:
                    fw.write(numpy.uint8(1).tobytes())
        if addSharp:
            fw.write(numpy.uint8(0).tobytes())

poly = scipy.interpolate.lagrange(numpy.array([78,65,67,71,84]), numpy.arange(2,7))

numerize_concat_fasta("/home/ljw/hg19_with_bowtie2_index/test.fa", poly, True, True, True, True)

for i in [0,2,4,5]:
    numerize_concat_fasta(f"/home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/local_file{i}", poly, False, False, False, False)
for i in [1,3,6,7]:
    numerize_concat_fasta(f"/home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/global_file{i}", poly, True, False, True, True)