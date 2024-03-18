#!/bin/bash

while IFS=, read -r path readlen cov; do
	b=$(basename ${path})
	s=${b%%.*}
	frac=$(awk -vcov=${cov} 'BEGIN { print 1 / cov }')
	echo "${s}"
	if [ ! -f coverage/${s}_cov.bed ]; then
		(samtools view -hb --region-file hg38_chroms.bed -s ${frac} --subsample-seed 123 ${path} | samtools depth - | python3 windows.py ${s}) &
	else
		echo "${s} already done"
	fi
	
	jobs=($(jobs -p))
	while (( ${#jobs[*]} >= 8 ))
	do
		sleep 30
		jobs=($(jobs -p))
	done
done < sample_list.csv
