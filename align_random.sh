#!/bin/bash
excu="/home/ljw/new_fold/old_desktop/shoujia/Graph_Projects/GAP/build/GAP"


threads_sz=24
max_extract=4;
diff_thres=1;
block_size=5
max_range=5;
max_mega=20000;




cd "/home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2"
mr=5

for round in $(seq 0 $mr); do
	rm -r tmp$round
done

Ti=-5
dT=-1
for round in $(seq 0 $mr); do
	mkdir tmp$round
	cd tmp$round
	if [ $round == 0 ]; then
		read_file="../read_file"
	else
		read_file="../tmp$(($round - 1))/fail"
	fi
	echo $read_file

	if [ ! -s $read_file ]; then
		break
	fi
	T=$(($Ti + $dT * $round))
	$excu --threads_sz $threads_sz --input $read_file --max_extract $max_extract --diff_thres $diff_thres --block_size $block_size --max_range $max_range --max_mega $max_mega --nodes node0,1,0 node1,1,0 node2,0,0 node3,0,0 node4,0,1 --shorts /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/local_file0.False.False.False.False.concat,node1,node1,,,,,$T /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/local_file2.False.False.False.False.concat,node1,node3,,,,,$T /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/local_file4.False.False.False.False.concat,node3,node1,,,,,$T /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/local_file5.False.False.False.False.concat,node2,node3,,,,,$T --longs /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/global_file1.True.False.True.True.concat,node0,node0,,,,,$T /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/global_file3.True.False.True.True.concat,node2,node3,,,,,$T /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/global_file6.True.False.True.True.concat,node3,node4,,,,,$T /home/ljw/new_fold/old_desktop/shoujia/test_RandomReads2/global_file7.True.False.True.True.concat,node0,node2,,,,,$T
	cd ..
done