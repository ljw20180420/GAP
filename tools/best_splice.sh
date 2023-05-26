#!/bin/bash

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	echo "choose the best indel results from all alternative inputs"
	echo "bash best_splice.sh output alternative_indel1 alternative_indel2 alternative_indel3 .."
	echo "is_better: prefer indels with at least one blunt end, and then choose the indel with the minimal indel (regardless of deletion and insertion)"
	echo "is_better2: the same as is_better but not prefer blunt end"
	echo "is_better3: only concern cut2 and key2, and then prefer the indels without insersion at cut2, and then choose the indel with the minimal deletion"
	exit 0
fi

output=$1;
shift 1;
fileseq=$@;
awk -v fileseq="$fileseq" '
function is_better(indels1, indels2){
	if (indels1[16]=="*")
		return 0;
	else
	{
		if (indels2[16]=="*")
			return 1;
		else
		{
			if (indels2[16]+indels2[18]>0 && indels2[17]+indels2[19]>0)
			{
				if (indels1[16]+indels1[18]==0 || indels1[17]+indels1[19]==0)
					return 1;
				else
				{
					if (indels1[16]+indels1[18]+indels1[17]+indels1[19] < indels2[16]+indels2[18]+indels2[17]+indels2[19])
						return 1;
					else
						return 0;
				}
			}
			else
			{
				if (indels1[16]+indels1[18]>0 && indels1[17]+indels1[19]>0)
					return 0;
				else
				{
					if (indels1[16]+indels1[18]+indels1[17]+indels1[19] < indels2[16]+indels2[18]+indels2[17]+indels2[19])
						return 1;
					else
						return 0;
				}
			}
		}
	}
	
}
function is_better2(indels1, indels2){
	if (indels1[16]=="*")
		return 0;
	else
	{
		if (indels2[16]=="*")
			return 1;
		else
		{
			if (indels1[16]+indels1[18]+indels1[17]+indels1[19] < indels2[16]+indels2[18]+indels2[17]+indels2[19])
				return 1;
			else
				return 0;
		}
	}
	
}
function is_better3(indels1, indels2){
	if (indels1[16]=="*")
		return 0;
	else
	{
		if (indels2[16]=="*")
			return 1;
		else
		{
			if (indels1[19] < indels2[19])
				return 1;
			else
			{
				if (indels1[19] > indels2[19])
						return 0;
				else
				{
					if (indels1[18] < indels2[18])
						return 1;
					else
						return 0;
				}
			}
		}
	}
}
BEGIN{
	filenum = split(fileseq, files, " ");
}
{
	if ($1!=QNAME)
	{	
		if (NR > 1)
			print line_best;
		QNAME = $1;
		start_flag = 1;
	}
	for (i=16; i<=19; ++i)
		indels[i]=$i;
	if (start_flag > 0 || is_better3(indels, indels_best)>0)
	{
		line_best=$0;
		for (i=16; i<=19; ++i)
			indels_best[i]=indels[i];
	}
	for (i=2; i<=filenum; ++i)
	{
		getline < files[i];
		for (j=16; j<=19; ++j)
			indels[j]=$j;
		if (is_better3(indels, indels_best)>0)
		{
			line_best=$0;
			for (j=16; j<=19; ++j)
				indels_best[j]=indels[j];
		}
	}
	start_flag = 0;
}
END{
	print $line_best;
}
' $1 > $output
