#!/bin/bash

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo "analyze the chimeric alignments based on the given rearrangement assumption"
	echo "bash get_indel.sh align_dir output cut1chr cut1pos cut1strand cut2chr cut2pos cut2strand thres extension align_mode"
	echo "thres: the longest allowed distance between cuts and key ends (the shifted ends from cuts by deletions or templated insertions)"
	echo "extension: the two-side extension around cuts to get references for alignments"
	echo "align_mode: try to minimize the distance between cut1 and key1 (up), cut2 and key2 (down), or (cut1+cut2)/2 and (key1+key2)/2 (others)"
	exit 0
fi

dir=`echo "${1%/}"`;
rm -f $dir/$2".indel";
touch $dir/$2".indel";
seqfile=`find $dir -name "*.fa"`;
countfile=${seqfile:0:-2}"txt";
files=`find $dir -name "*.ext"`;
bedtools getfasta -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed <(printf "%s\t%d\t%d\t%s\t%s\t%s\n" $3 $(($4 - ${10})) $(($4 + ${10})) "cut1ref" "." $5) | tail -n1 > cut1ref;
bedtools getfasta -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed <(printf "%s\t%d\t%d\t%s\t%s\t%s\n" $6 $(($7 - ${10})) $(($7 + ${10})) "cut1ref" "." $8) | tail -n1 > cut2ref;
for file in $files; do
	awk -v alg=${file:0:-3}"alg" -v countfile=$countfile -v seqfile=$seqfile -v cut1chr=$3 -v cut1pos=$4 -v cut1strand=$5 -v cut2chr=$6 -v cut2pos=$7 -v cut2strand=$8 -v thres=$9 -v cutrefpos=${10} -v align_mode=${11} '
	function shift_limit(cut1refarray, cut1reflen, cut2refarray, cut2reflen, map1refkey, map2refkey, downup)
	{
		j=0;
		for(i=0; i<map1refkey && j<map2refkey; ++i)
		{
			if (cut1refarray[map1refkey-i]!=cut2refarray[map2refkey-j])
				break;
			++j;
		}
		downup[1]=-j;
		j=1;
		for(i=1; i<=cut1reflen-map1refkey && j<=cut2reflen-map2refkey; ++i)
		{
			if (cut1refarray[map1refkey+i]!=cut2refarray[map2refkey+j])
				break;
			++j;
		}
		downup[2]=j-1;
	}
	function get_key(isfirst, strand, mapstart, mapend, refstart, refend, mapkey){
		if (isfirst>0)
		{
			if (strand=="+")
			{
				mapkey[1] = mapend;
				maprefkey = mapkey[1] - refstart;
			}
			else
			{
				mapkey[1] = mapstart;
				maprefkey = refend - mapkey[1];
			}
		}
		else
		{
			if (strand=="+")
			{
				mapkey[2] = mapstart;
				maprefkey = mapkey[2] - refstart;
			}
			else
			{
				mapkey[2] = mapend;
				maprefkey = refend - mapkey[2];
			}
		}
		return maprefkey;
	}
	function get_indel(map1refkey, map2refkey, downup, cut1refpos, cut2refpos, queryfull, query1end, query2start, indels, align_mode){
		indels[5]=substr(queryfull, query1end+1, query2start-query1end);
		if (align_mode == "up")
			target_shift = cut1refpos - map1refkey;
		else
		{
			if (align_mode == "down")
				target_shift = cut2refpos - map2refkey;
			else
				target_shift = int((cut1refpos + cut2refpos - map1refkey - map2refkey)/2);
		}
		if (target_shift < downup[1])
			target_shift = downup[1];
		if (target_shift > downup[2])
			target_shift = downup[2];
		map1refkey += target_shift;
		map2refkey += target_shift;
		if (map1refkey < cut1refpos)
		{
			indels[1] = cut1refpos - map1refkey;
			indels[2] = 0;
		}
		else
		{
			indels[1] = 0;
			indels[2] = map1refkey - cut1refpos;
		}
		if (map2refkey < cut2refpos)
		{
			indels[3] = 0;
			indels[4] = cut2refpos - map2refkey; 
		}
		else
		{
			indels[3] = map2refkey - cut2refpos;
			indels[4] = 0;
		}
		return target_shift;
	}
	function getRC(seq, RCmap){
		n = split(seq, seqarray, "");
		seqRC="";
		for (i=1; i<=n; ++i)
			seqRC = RCmap[seqarray[i]] seqRC;
		return seqRC;
	}
	BEGIN{
		getline cut1ref < "cut1ref";
		getline cut2ref < "cut2ref";
		cut1ref = toupper(cut1ref);
		cut2ref = toupper(cut2ref);
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
		
		n = split(cut1ref, cut1refarray, "");
		n = split(cut2ref, cut2refarray, "");

		getline name < seqfile;
		if (name ~ /^@/)
			mode = "fastq";
		else
			mode = "fasta";
		getline queryfull < seqfile;
		if (mode == "fastq")
		{
			getline junk < seqfile;
			getline junk < seqfile;
		}
		getline countline < countfile;
		split(countline, counts, " ");
		count = counts[3];
	}
	{
		getline < alg;
		getline query < alg;
                getline mid < alg;
		getline ref < alg;
		segnum=split(query,query_segs," ");
		split(ref,ref_segs," ");		

		while (name != $1)
		{
			getline name < seqfile;
			getline queryfull < seqfile;
			if (mode == "fastq")
			{
				getline junk < seqfile;
				getline junk < seqfile;
			}
			getline countline < countfile;
			split(countline, counts, " ");
			count = counts[3];
		}
		
		QNAME=substr($1, 2, length($1)-1);
		MAPIDX=$2;
		for (i=1; i<=segnum; ++i)
		{
			getline;
			n=split($1,a,":");
			if (a[n] ~ /_RC$/)
			{
				RNAMEs[i]=substr(a[n], 2, length(a[n])-4);
				strand[i]="-";
				mapend[i]=ref2len[RNAMEs[i]]-$2;
			}
			else
			{
				RNAMEs[i]=substr(a[n], 2);
				strand[i]="+";
				mapstart[i]=$2;
			}
			query_start[i]=$(NF-1);
			score[i]=$NF;
			getline;
			if (a[n] ~ /_RC$/)
				mapstart[i]=ref2len[RNAMEs[i]]-$2;
			else
				mapend[i]=$2;
			query_end[i]=$(NF-1);
			score[i]=$NF-score[i];
		}

		map1refkey = get_key(1, strand[1], mapstart[1], mapend[1], cut1pos-cutrefpos, cut1pos+cutrefpos, mapkey);
		if (segnum > 1)
			map2refkey = get_key(0, strand[segnum], mapstart[segnum], mapend[segnum], cut2pos-cutrefpos, cut2pos+cutrefpos, mapkey);
		else
			mapkey[2]="*";

		if (segnum<2 || RNAMEs[1]!=cut1chr || RNAMEs[segnum]!=cut2chr || strand[1]!=cut1strand || strand[segnum]!=cut2strand)
		{
			for (i=1; i<=5; ++i)
				indels[i]="*";
			downup[1]="*";
			downup[2]="*";
			target_shift="*";
		}
		else
		{
			if (map1refkey-cutrefpos-thres>0 || map1refkey-cutrefpos+thres<0 || map2refkey-cutrefpos-thres>0 || map2refkey-cutrefpos+thres<0)
			{
				for (i=1; i<=5; ++i)
					indels[i]="*";
				downup[1]="*";
				downup[2]="*";
				target_shift="*";
			}
			else
			{
				if (query_start[segnum] > query_end[1])
				{
					downup[1]=0;
					downup[2]=0;
				}
				else
					shift_limit(cut1refarray, length(cut1ref), cut2refarray, length(cut2ref), map1refkey, map2refkey, downup);
				target_shift=get_indel(map1refkey, map2refkey, downup, cutrefpos, cutrefpos, queryfull, query_end[1], query_start[segnum], indels, align_mode);
			}
		}
		if (segnum<2)
		{
			strand_sec="*";
			RNAME_sec="*"
			ref_seg_sec="*";
			query_seg_sec="*";
			mapstart_sec="*";
			mapend_sec="*";		
		}
		else
		{
			strand_sec=strand[segnum]; 
			RNAME_sec=RNAMEs[segnum];
			ref_seg_sec=ref_segs[segnum];
			query_seg_sec=query_segs[segnum];
			mapstart_sec=mapstart[segnum];
			mapend_sec=mapend[segnum];	
		}
		printf("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", QNAME, count, MAPIDX, RNAMEs[1], mapkey[1], strand[1], RNAME_sec, mapkey[2], strand_sec, cut1chr, cut1pos, cut1strand, cut2chr, cut2pos, cut2strand, indels[1], indels[2], indels[3], indels[4], indels[5], downup[1], downup[2], target_shift, ref_segs[1], query_segs[1], mapstart[1], mapend[1], ref_seg_sec, query_seg_sec, mapstart_sec, mapend_sec);
		
	}
	' $file >> $dir/$2".indel";
done
rm -f cut1ref cut2ref
