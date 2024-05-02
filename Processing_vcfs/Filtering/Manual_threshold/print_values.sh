#!/bin/bash

### script for returning vcf lines for directory of gw plotted SVs

vcfdir="/path/to/project/variants/hawk/filtered"
flaseneg="/path/to/project/variants/processing_vcfs/manual_filtering/hawk/DEL/FN"
for fn in ${flaseneg}/*.png; do
	IFS="~" read -r -a loci <<< "${fn}"
	echo "${loci[@]}"
	for i in ${vcfdir}/*.vcf; do
		bcftools filter -i"CHROM == \"${loci[1]}\" && POS == $((${loci[2]} + 1))" ${i} > tmp_var.vcf
		if [[ $(tail -n1 tmp_var.vcf) != "#"* ]]; then
			echo "$(tail -n1 tmp_var.vcf)"
		fi
	done
done
