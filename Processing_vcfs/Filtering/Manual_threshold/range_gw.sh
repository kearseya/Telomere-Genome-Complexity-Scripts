#!/bin/bash


i=$(ls /path/to/project/variants/icgc/normal/$1*)
b=$(basename $i .vcf)

# bcftools filter -i"INFO/SVTYPE == \"$2\" && FORMAT/PROB <= $4 && FORMAT/PROB >= $3" /path/to/project/variants/icgc/normal/${b}.vcf | tail -n1

echo "$2 $3-$4"
bcftools filter -i"INFO/SVTYPE == \"$2\" && FORMAT/PROB <= $4 && FORMAT/PROB >= $3" /path/to/project/variants/icgc/normal/${b}.vcf | bcftools sort - | ~/Applications_Downloaded/tmpforks/gw/./gw --link sv --labels yes,no,maybe /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/${b}.bam -v - --track ${i}
