#!/bin/bash

### script for generating SVs images, putting them into folders for manual classification into true/false positive/negative
### either with large thumbnail or gnome-sushi/preview to checkout SVs with ctrl-click dragging into correct folders
### use range_gw.sh first to discover ranges for suitable probability values

vcfdir="/path/to/project/variants/hawk/filtered"
types=("DEL" "DUP" "INV" "TRA" "INS")
declare -A probs
probs=( ["DEL"]=0.15 ["DUP"]=0.2 ["INV"]=0.45 ["TRA"]=0.45 ["INS"]=0.45 )
outdir="/path/to/project/variants/processing_vcfs/manual_filtering/hawk"


## set up directory structure
subdirs=("TP" "TN" "FP" "FN")
for t in ${types[@]}; do
	mkdir -p "${outdir}/${t}/${probs[${t}]}/above/rand"
	mkdir -p "${outdir}/${t}/${probs[${t}]}/below/rand"
	for s in ${subdirs[@]}; do
		mkdir -p "${outdir}/${t}/${probs[${t}]}/${s}"
	done
done

## generate sv images with gw
for i in ${vcfdir}/*.vcf; do
	s=$(basename "${i}" .vcf)
	b=$(ls /mnt/breast/${s}.*am)
	echo "vcf ${i}\nbase ${s}\nbam ${b}"
	for t in ${types[@]}; do
		if [ ${t} != "INS" ] && [ ${t} != "TRA" ]; then
			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB >= ${probs[${t}]} && INFO/SVLEN >= 10000" ${i} > tmp.vcf
			gw --bam "${b}" -v tmp.vcf -o "${outdir}/${t}/${probs[${t}]}/above" -n /scratch/ProjectDir/hg38.fa
			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB < ${probs[${t}]} && INFO/SVLEN >= 10000" ${i} > tmp.vcf
			gw --bam "${b}" -v tmp.vcf -o "${outdir}/${t}/${probs[${t}]}/below" -n /scratch/ProjectDir/hg38.fa
		else
			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB >= ${probs[${t}]}" ${i} > tmp.vcf
			gw --bam "${b}" -v tmp.vcf -o "${outdir}/${t}/${probs[${t}]}/above" -n /scratch/ProjectDir/hg38.fa
			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB < ${probs[${t}]}" ${i} > tmp.vcf
			gw --bam "${b}" -v tmp.vcf -o "${outdir}/${t}/${probs[${t}]}/below" -n /scratch/ProjectDir/hg38.fa
		fi
	done
done

## create random subsamples
random=true
categories=("above" "below")
if [ ${random} == true ]; then
	N=50 # number from each category
	for t in ${types[@]}; do
		for a in ${categories[@]}; do
			ls "${outdir}/${t}/${probs[${t}]}/${a}/" | sort -R | tail -$N | while read file; do
				mv "${outdir}/${t}/${probs[${t}]}/${a}/${file}" "${outdir}/${t}/${probs[${t}]}/${a}/rand/${file}"
			done
		done
	done
fi


### streaming doesnt work for making images
# for i in ${vcfdir}/*.vcf; do
# 	s=$(basename "${i}" .vcf)
# 	b=$(ls /mnt/breast/${s}.*am)
# 	echo "vcf ${i}\nbase ${s}\nbam ${b}"
# 	for t in ${types[@]}; do
# 		if [ ${t} != "INS" ]; then
# 			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB >= ${probs[${t}]} && INFO/SVLEN >= 10000" ${i} | gw --bam "${b}" -v - -o "${outdir}/${t}/mt${probs[${t}]}" /scratch/ProjectDir/hg38.fa
# 			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB < ${probs[${t}]} && INFO/SVLEN >= 10000" ${i} | gw --bam "${b}" -v - -o "${outdir}/${t}/lt${probs[${t}]}" /scratch/ProjectDir/hg38.fa
# 		else
# 			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB >= ${probs[${t}]}" ${i} | gw --bam "${b}" -v - -o "${outdir}/${t}/mt${probs[${t}]}" -n /scratch/ProjectDir/hg38.fa
# 			bcftools filter -i"INFO/SVTYPE == \"${t}\" && FORMAT/PROB < ${probs[${t}]}" ${i} | gw --bam "${b}" -v - -o "${outdir}/${t}/lt${probs[${t}]}" -n /scratch/ProjectDir/hg38.fa
# 		fi
# 	done
# done


