#!/bin/bash

module load samtools 

in_dir=/scratch/ProjectDir/aligns
ref=/scratch/ProjectDir/reference_genomes/hg38.fa

for i in ${in_dir}/*.cram;
do
    s=$(basename ${i} .cram);
    echo ${s};
    samtools depth -a -r chr22 ${i} > ${s}.chr22.tsv;
    /scratch/User/miniconda3/bin/python3 gc_tabulate.py -d ${s}.chr22.tsv -ref ${ref} > ${s}_gc_cov.tsv;
    rm ${s}.chr22.tsv;
done

