#!/bin/bash

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo "filter the same file according to chromosome, range, orientation, and gap"
	echo "only those alignments with exactly two segments are remained"
	echo "the first segment must be upstream of the second one"
	echo "bash filter_sam.sh dir chromosome locus1 locus2 orientation gap"
	echo "chromosome: both segments map to this chromosome"
	echo "locus1<locus2 define the view window"
	echo "orientation is of the form ** (the first and second strands), where *=a,+,-. a means the strand does not matter"
	echo "gap is the minimal endurable distance (gap) between the first and second segments"
	exit 0
fi

dir=`echo "${1%/}"`;
sam=`find $dir -name "*.sam" | grep -v "filter" | grep -v "."$4"."$5"."`;
awk -v chr=$2 -v locus1=$3 -v locus2=$4 -v ori=$5 -v gap=$6 '
function CIGAR_len(CIGAR){
	n = split(CIGAR, array, /[M=XID]/, seps);
	len=0;
	for (i=1;i<n;++i)
		if(seps[i]!="D")
			len+=array[i];
	return len;
}
{
	if(QNAME!=$1)
	{
		if (ori ~ /^a/)
			ori_left=1;
		else
		{
			ori_left=0;
			if (ori ~ /^+/ && int(flags[1]/16)%2==0 || ori ~ /^-/ && int(flags[1]/16)%2==1)
				ori_left=1;
		}
		if (ori ~ /a$/)
			ori_right=1;
		else
		{
			ori_right=0;
			if (segnum > 1 && (ori ~ /+$/ && int(flags[2]/16)%2==0 || ori ~ /-$/ && int(flags[2]/16)%2==1))
				ori_right=1;
		}
		if (segnum==2 && RNAMEs[1]==chr && RNAMEs[2]==chr && ori_left>0 && ori_right>0 && first_right<=second_left && first_right > locus1 && second_left <= locus2 && second_left - first_right > gap)
			for (i = 1; i <= segnum; ++i)
				print lines[i];
		QNAME=$1;
		segnum=0;
	}
	++segnum;
	lines[segnum]=$0;
	flags[segnum]=$2;
	RNAMEs[segnum]=$3;
	if (segnum==1)
		first_right = $4 + CIGAR_len($6);
	else if (segnum==2)
		second_left =$4;
}
END{
	if (segnum==2 && RNAMEs[1]==chr && RNAMEs[2]==chr && ori_left>0 && ori_right>0 && first_right<=second_left && first_right > locus1 && second_left <= locus2 && second_left - first_right > gap)
			for (i = 1; i <= segnum; ++i)
				print lines[i];
}
' $sam > ${sam:0:-3}filter.$5.$6.sam;
samtools view -b -T /home/ljw/hg19_with_bowtie2_index/hg19.fa ${sam:0:-3}filter.$5.$6.sam | samtools sort -o ${sam:0:-3}filter.$5.$6.bam;
samtools index ${sam:0:-3}filter.$5.$6.bam;
