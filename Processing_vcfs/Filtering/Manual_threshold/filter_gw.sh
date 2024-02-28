#!/bin/bash

# bcftools filter -e'(SVTYPE == "DEL" && PROB > 0.7) || (SVTYPE == "INS" && PROB > 0.7) || (SVTYPE == "TRA" && PROB > 0.3) || (SVTYPE == "DUP" && PROB > 0.1) || (SVTYPE == "INV" && PROB > 0.25)' /path/to/project/variants/icgc/normal/02226b857f74317cd2cd8a389b186e0c.vcf | gw --link sv --labels yes,no,maybe --out-lables all_types.tsv /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/02226b857f74317cd2cd8a389b186e0c.bam -v -

i=$(ls /path/to/project/variants/icgc/normal/$1*)
b=$(basename $i .vcf)

echo "Deletions 0.20 - 0.40"
bcftools filter -i'INFO/SVTYPE == "DEL" && FORMAT/PROB <= 0.55 && FORMAT/PROB <= 0.3' /path/to/project/variants/icgc/normal/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/${b}.bam -v - --track ${i}

echo "Insertions 0.0 - 0.20"
bcftools filter -i'INFO/SVTYPE == "INS" && FORMAT/PROB <= 0.2 && FORMAT/PROB >= 0' /path/to/project/variants/icgc/normal/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/${b}.bam -v - --track ${i}

echo "Duplications 0.0 - 0.25"
bcftools filter -i'INFO/SVTYPE == "DUP" && FORMAT/PROB <= 0.35 && FORMAT/PROB >= 0.1' /path/to/project/variants/icgc/normal/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/${b}.bam -v - --track ${i}

echo "Inversions 0.0 - 0.25"
bcftools filter -i'INFO/SVTYPE == "INV" && FORMAT/PROB <= 0.35 && FORMAT/PROB <= 0.1' /path/to/project/variants/icgc/normal/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/${b}.bam -v - --track ${i}

echo "Translocations 0.0 - 0.35"
bcftools filter -i'INFO/SVTYPE == "TRA" && FORMAT/PROB <= 0.45 && FORMAT/PROB <= 0.1' /path/to/project/variants/icgc/normal/${b}.vcf | gw --link sv --labels yes,no,maybe --out-labels $1_del.tsv /path/to/project/TL_prediction/reference_genomes/references_hs37d5_hs37d5.fa -b /mnt/icgc/trimmed_bam/${b}.bam -v - --track ${i}
