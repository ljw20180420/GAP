#!/bin/bash

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo "transform alg to sam"
	echo "for matches, = but not M is used"
	echo "bash to_sam.sh input.alg.dir output.sam"
	exit 0
fi 

dir=`echo "${1%/}"`;
rm -f $dir/$2".sam";
touch $dir/$2".sam";
files=`find $dir -name "*.ext"`;
for file in $files; do
	awk -v alg=${file:0:-3}"alg" '
	function getCIGARs(QNAME, query_segs, mid, ref_segs, segnum, CIGARs){
		start=1;
                for (i=1; i<=segnum; ++i)
		{
			ql = length(query_segs[i])
			mid_segs[i]=substr(mid,start,length(query_segs[i]));
			sub(/^-+/, "", query_segs[i]);
			mid_segs[i]=substr(mid_segs[i], length(mid_segs[i])-length(query_segs[i])+1, length(query_segs[i]));
			ref_segs[i]=substr(ref_segs[i], length(ref_segs[i])-length(query_segs[i])+1, length(query_segs[i]));
			sub(/-+$/, "", query_segs[i]);
			mid_segs[i]=substr(mid_segs[i], 1, length(query_segs[i]));
			ref_segs[i]=substr(ref_segs[i], 1, length(query_segs[i]));
			m = split(query_segs[i], query_array, "");
			split(mid_segs[i], mid_array, "");
			split(ref_segs[i], ref_array, "");
			typepre="?";
			cumsum=0;
			CIGARs[i]="";
			for (j=1; j<=m; ++j)
			{
				if (mid_array[j]=="|")
					type="=";
				else
				{
					if (query_array[j]=="-")
						type="D";
					else
					{
						if (ref_array[j]=="-")
							type="I";
						else
							type="X";
					}
				}
				if (type != typepre)
				{
					if (cumsum > 0)
					{
						CIGARs[i] = CIGARs[i] cumsum typepre; 
						cumsum = 0;
					}
					typepre=type;
				}
				++cumsum;
			}
			CIGARs[i] = CIGARs[i] cumsum type;
			start += ql+1;
		}
	}
	function getFLAGs(segnum, FLAGs, strand){
		for (i=1; i<=segnum; ++i)
		{
			FLAGs[i] = 2;
			if (segnum > 1)
				FLAGs[i] += 1;
			if (strand[i]=="-")
				FLAGs[i] += 16;
			if (strand[i%segnum + 1]=="-")
				FLAGs[i] += 32;
			if (i!=1)
				FLAGs[i] += 2048;
		}
	}
	function getTAGs(segnum, TAGs, RNAMEs, ref_start, strand, CIGARs){
		for (i=1; i<=segnum; ++i)
		{
			TAGs[i] = sprintf("XX:i:%d\tSA:Z:", i);
			for (j=1; j<=segnum; ++j)
				if (j!=i)
					TAGs[i] = sprintf("%s%s,%s,%s,%s,*,0;", TAGs[i], RNAMEs[j], ref_start[j], strand[j], CIGARs[j]);
		}
	}
	function getTLENs(segnum, TLENs, RNAMEs, ref_start, ref_end)
	{
		if (segnum == 1)
			TLENs[1]=0;
		else
		{
			for (i=2; i<=segnum; ++i)
				if (RNAMEs[i]!=RNAMEs[1])
				{
					for (j=1; j<=segnum; ++j)
						TLENs[j]=0;
					return;
				}
			leftmost = ref_start[1];
			rightmost = ref_end[1];
			rightidx = 1;
			for (i=2; i<=segnum; ++i)
			{
				if (leftmost > ref_start[i])
					leftmost = ref_start[i];
				if (rightmost < ref_end[i])
				{
					rightmost = ref_end[i];
					rightidx = i;
				}
			}
			TLEN = rightmost - leftmost;
			for (i=1; i<=segnum; ++i)
				TLENs[i] = ""TLEN;
			TLENs[rightidx] = "-"TLENs[rightidx];
		}
	}
	BEGIN{
		while((getline < "/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes")>0)
			ref2len[$1]=$2;
		RCmap["A"]="T";
		RCmap["a"]="T";
		RCmap["T"]="A";
		RCmap["t"]="A";
		RCmap["G"]="C";
		RCmap["g"]="C";
		RCmap["C"]="G";
		RCmap["c"]="G";
		RCmap["N"]="N";
		RCmap["n"]="N";
	}
	{
		getline < alg;
		getline query < alg;
                getline mid < alg;
		getline ref < alg;
		segnum=split(query,query_segs," ");
		split(ref,ref_segs," ");

		QNAME=substr($1, 2, length($1)-1);
		for (i=1; i<=segnum; ++i)
		{
			getline;
			n=split($1,a,":");
			if (a[n] ~ /_RC$/)
			{
				RNAMEs[i]=substr(a[n], 2, length(a[n])-4);
				strand[i]="-";
				ref_end[i]=ref2len[RNAMEs[i]]-$2+1;
			}
			else
			{
				RNAMEs[i]=substr(a[n], 2);
				strand[i]="+";
				ref_start[i]=$2+1;
			}
			query_start[i]=$(NF-1);
			score[i]=$NF;
			getline;
			if (a[n] ~ /_RC$/)
				ref_start[i]=ref2len[RNAMEs[i]]-$2+1;
			else
				ref_end[i]=$2+1;
			query_end[i]=$(NF-1);
			score[i]=$NF-score[i];
		}
		
		getCIGARs(QNAME, query_segs, mid, ref_segs, segnum, CIGARs);
		getFLAGs(segnum, FLAGs, strand);
		getTAGs(segnum, TAGs, RNAMEs, ref_start, strand, CIGARs);
		getTLENs(segnum, TLENs, RNAMEs, ref_start, ref_end);
		for (i=1; i<=segnum; ++i)
		{
			gsub(/-/, "", query_segs[i]);
			if (strand[i]=="+")
				SEQ=query_segs[i];
			else
			{
				m=split(query_segs[i],b,"");
				SEQ="";
				for (j=1; j<=m; ++j)
					SEQ = RCmap[b[j]] SEQ;
			}
			printf("%s\t%d\t%s\t%d\t255\t%s\t%s\t%d\t%s\t%s\t*\t%s\n", QNAME, FLAGs[i], RNAMEs[i], ref_start[i], CIGARs[i], RNAMEs[i%segnum + 1], ref_start[i%segnum + 1], TLENs[i], SEQ, TAGs[i]);
		}
	}
	' $file >> $dir/$2".sam";
done
samtools view -b -T /home/ljw/hg19_with_bowtie2_index/hg19.fa $dir/$2".sam" | samtools sort -o $dir/$2".bam";
samtools index $dir/$2".bam";
