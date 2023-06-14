#!/bin/bash

excu="/home/ljw/new_fold/old_desktop/shoujia/Graph_Projects/GAP/build/GAP"
genome="/home/ljw/hg19_with_bowtie2_index/hg19.fa"
msn=1
bs=5
tz=24
cd "/home/ljw/new_fold/old_desktop/shoujia/test_time_loop3"
mr=5

for round in $(seq 0 $mr); do
	rm -r tmp$round
done
rm -r tmplast

Ti=-5
dT=-5
for round in $(seq 0 $mr); do
	mkdir tmp$round
	cd tmp$round
	if [ $round == 0 ]; then
		read_file="../read_file.fa"
	else
		read_file="../tmp$(($round - 1))/fail"
	fi
	echo $read_file

	if [ ! -s $read_file ]; then
		break
	fi
	T=$(($Ti + $dT * $round))
	$excu ---threads_sz $tz ---read $read_file ---min_seg_num 1 ---max_seg_num $msn ---block_size $bs ---nodes --names node0 ---roots node0 ---targets node0 ---globals --names $genome --tails node0 --heads node0 --Ts $T --r "true"
	cd ..
done
mkdir tmplast
cd tmplast
read_file="../tmp$mr/fail"
echo $read_file
if [ -s $read_file ]; then
	$excu ---threads_sz $tz ---read $read_file ---block_size $bs ---nodes --names $(seq -f 'node%.f' -s ' ' 0 $msn) ---roots node0 ---targets node$msn ---globals --names $(yes $genome | head -n $msn | tr "\n" " ") --tails $(seq -f 'node%.f' -s ' ' 0 $(($msn - 1))) --heads $(seq -f 'node%.f' -s ' ' 1 $msn) --Ts 0 --r $(yes true | head -n $msn | tr "\n" " ")
fi

