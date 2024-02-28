#!/bin/bash

while read l; do
	echo $l
	#((p=p+1))
	p=$(echo $l | cut -d, -f1)
	t=$(echo $l | cut -d, -f3)
	n=$(echo $l | cut -d, -f2)
	echo "PID: ${p}, Normal: ${n}, Tumour: ${t}"
	telomerehunter -ibt /mnt/breast/${t}.*am -ibc /mnt/breast/${n}.*am -o telomerehunter_workdir -p ${p} -b /path/to/project/variants/processing_vcfs/hg38_cytoBand.txt -pl
done < sample_pairs_index.csv

