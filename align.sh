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
	$excu --threads_sz $tz --input $read_file --min_seg_num 1 --max_seg_num $msn --block_size $bs --nodes node0,1,1 --longs $genome,node0,node0,,,,,$T
	cd ..
done
mkdir tmplast
cd tmplast
read_file="../tmp$mr/fail"
echo $read_file

nodes=""
for i in $(seq 0 $msn); do
	if [ $i == 0 ]; then
		is_root=1
	else
		is_root=0
	fi
	if [ $i == $msn ]; then
		is_target=1
	else
		is_target=0
	fi
	nodes=$nodes" node"$i","$is_root","$is_target
done
echo $nodes

longs=""
for i in $(seq 0 $(($msn-1))); do
	longs=$longs" "$genome",node"$i",node"$(($i+1))",,,,,0.0"
done
echo $longs

if [ -s $read_file ]; then
	$excu --threads_sz $tz --input $read_file --block_size $bs --nodes $nodes --longs $longs
fi