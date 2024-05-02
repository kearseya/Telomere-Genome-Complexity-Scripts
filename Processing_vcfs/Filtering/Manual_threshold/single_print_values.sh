#!/bin/bash

### script to collect vcf line from gw SV plot, drag png into terminal for filepath

vcfdir="/path/to/project/variants/hawk/filtered"
fn=$(basename $1 .png)
#loci=(${${fn}//~/ })
#flaseneg="/path/to/project/variants/processing_vcfs/manual_filtering/hawk/DEL/FN"
IFS="~" read -r -a loci <<< "${fn}"
echo "${loci[@]}"
for i in ${vcfdir}/*.vcf; do
	bcftools filter -i"CHROM == \"${loci[1]}\" && POS == $((${loci[2]} + 1))" ${i} > tmp_var.vcf
	if [[ $(tail -n1 tmp_var.vcf) != "#"* ]]; then
		echo "$(tail -n1 tmp_var.vcf)"
		#s=$(basename ${i} .vcf)
		#b=$(ls /mnt/breast/${s}.*am)
		#gw --bam "${b}" -v tmp_var.vcf /scratch/ProjectDir/hg38.fa
	fi
done
