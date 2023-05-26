#!/bin/bash

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo "count and classify the reads in fastq"
	echo "bash percent.sh input.fastq output percentage_threshold"
	exit 0
fi

total=$(wc -l $1 | cut -f1 -d" ");
awk '{if(NR%4 == 2) print $0;}' $1 | sort -f | uniq -ic | sort -bnr -k1,1 | awk -v tot=$[total/4] -v thres=$3 '{per=$1/tot; if(per > thres) print $2,per,$1;}' > $2;
