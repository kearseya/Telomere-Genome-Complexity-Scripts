#!/bin/bash

mkdir -p list_dir
mkdir -p log_dir
mkdir -p jasmine_filtered
mkdir -p merged_dir

for i in /path/to/project/variants/bgi/tumour/*.vcf;
do
	b=$(basename ${i} .vcf)
	echo "Processing ${b}"
	cat normals.txt > list_dir/${b}.txt
	echo ${i} >> list_dir/${b}.txt
	jasmine --mark_specific file_list=list_dir/${b}.txt out_file=${b}_merge.vcf > log_dir/${b}.log
	mv output/${b}_markedSpec.vcf jasmine_filtered/
	mv ${b}_merge.vcf merged_dir
done
