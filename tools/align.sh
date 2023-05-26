#!/bin/bash

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo "align reads to the genome by genome-wide repeated alignment in 5 rounds (0~4) with increasing basic penalty T (default)"
	echo "the final step align reads by consecutive cross edges up to maximally allowed segments"
	echo "before running, copy GAP excutable to the running dir"
	echo "bash align.sh dir max_seg_num genome max_round T dT block_size excutable"
	exit 0
fi

dir=${1%/}
msn=${2:-"1"}
genome=${3:-"/home/ljw/hg19_with_bowtie2_index/hg19.fa"}
mr=${4:-"4"}
T=${5:-"-10"}
dT=${6:-"-5"}
bs=${7:-"100"}
excu=${8:-"GAP"}
tz=${9:-"24"}

run=$(echo ${dir##*/})
reads=""
for suffix in {fa,fq,fasta,fastq}; do
	reads=$reads" "$(find $dir -name "*.$suffix")
done
cd $dir;
$excu ---threads_sz $tz ---run $run ---reads $reads ---max_round $mr ---min_seg_num 1 ---max_seg_num $msn ---block_size $bs ---nodes --names node0 ---roots node0 ---targets node0 \
---globals --names /home/ljw/hg19_with_bowtie2_index/hg19.fa --tails node0 --heads node0 --Ts $T --dTs $dT --r "true";
if [ -s $run$mr".fail" ]; then
	$excu ---threads_sz $tz ---run ${run}_last ---reads $run$mr.fail ---block_size $bs ---nodes --names $(seq -f 'node%.f' -s ' ' 0 $msn) ---roots node0 ---targets node$msn \
	---globals --names $(yes $genome | head -n $msn | tr "\n" " ") \
	--tails $(seq -f 'node%.f' -s ' ' 0 $(($msn - 1))) --heads $(seq -f 'node%.f' -s ' ' 1 $msn) --Ts 0 --r $(yes true | head -n $msn | tr "\n" " ")
fi

