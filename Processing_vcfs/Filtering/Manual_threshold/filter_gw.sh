#!/bin/bash

vcfdir="/path/to/vcfs"
bamdir="/path/to/bams"
refpath="/path/to/reference.fa"

i=$(ls ${vcfdir}/$1*)
b=$(basename $i .vcf)

echo "Deletions 0.20 - 0.40"
bcftools filter -i'INFO/SVTYPE == "DEL" && FORMAT/PROB <= 0.55 && FORMAT/PROB <= 0.3' ${vcfdir}/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv ${refpath} -b ${bamdir}/${b}.bam -v - --track ${i}

echo "Insertions 0.0 - 0.20"
bcftools filter -i'INFO/SVTYPE == "INS" && FORMAT/PROB <= 0.2 && FORMAT/PROB >= 0' ${vcfdir}/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv ${refpath} -b ${bamdir}/${b}.bam -v - --track ${i}

echo "Duplications 0.0 - 0.25"
bcftools filter -i'INFO/SVTYPE == "DUP" && FORMAT/PROB <= 0.35 && FORMAT/PROB >= 0.1' ${vcfdir}/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv ${refpath} -b ${bamdir}/${b}.bam -v - --track ${i}

echo "Inversions 0.0 - 0.25"
bcftools filter -i'INFO/SVTYPE == "INV" && FORMAT/PROB <= 0.35 && FORMAT/PROB <= 0.1' ${vcfdir}/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv ${refpath} -b ${bamdir}/${b}.bam -v - --track ${i}

echo "Translocations 0.0 - 0.35"
bcftools filter -i'INFO/SVTYPE == "TRA" && FORMAT/PROB <= 0.45 && FORMAT/PROB <= 0.1' ${vcfdir}/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv ${refpath} -b ${bamdir}/${b}.bam -v - --track ${i}
