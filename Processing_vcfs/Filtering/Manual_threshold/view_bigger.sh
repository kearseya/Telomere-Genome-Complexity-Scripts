#!/bin/bash

### script to find and open SV in gw from a previously plotted png
### drag png into terminal for the filepath
### (used because sample names not in file names, plotted collective)

vcfdir="/path/to/project/variants/hawk/filtered"
fn=$(basename $1 .png)
IFS="~" read -r -a loci <<< "${fn}"
echo "${loci[@]}"

for i in ${vcfdir}/*.vcf; do
	bcftools filter -i"CHROM == \"${loci[1]}\" && POS == $((${loci[2]} + 1))" ${i} > tmp_var.vcf
	if [[ $(tail -n1 tmp_var.vcf) != "#"* ]]; then
		echo "$(tail -n1 tmp_var.vcf)"
		s=$(basename ${i} .vcf)
		b=$(ls /mnt/breast/${s}.*am)
		gw --bam "${b}" -v tmp_var.vcf /scratch/ProjectDir/hg38.fa
	fi
done
