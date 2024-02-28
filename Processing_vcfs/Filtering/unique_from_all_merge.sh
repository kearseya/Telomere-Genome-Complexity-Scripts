#!/bin/bash
normal=$1
tumour=$2

mkdir -p bash_filtered
mkdir -p bash_filtered/merged

for i in ${tumour}/*.vcf;
do
	b=$(basename ${i} .vcf)
	echo "merging ${b}"
	dysgu merge -o bash_filtered/merged/${b}_merged.vcf ${normal}/*.vcf ${i}
	echo "filtering ${b}"
	./unique_from_all_merge.py -o bash_filtered bash_filtered/merged/${b}_merged.vcf ${b}
done
