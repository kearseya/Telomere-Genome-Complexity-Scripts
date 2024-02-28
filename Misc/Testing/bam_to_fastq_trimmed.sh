#!/bin/bash

module load picard/2.20.2

n=216

java -jar /software/genomics/picard/2.20.2/INSTALL/picard.jar SamToFastq I=realigned_hg19/tmp/DB216_tel.bam FASTQ=fastq_trimmed/DB216_1.fq SECOND_END_FASTQ=fastq_trimmed/DB216_2.fq UNPAIRED_FASTQ=fastq_trimmed/DB216_3.fq

#for i in $(ls realigned_hg19/tmp/*bam | head); do n=$(basename ${i} .bam); java -jar /software/genomics/picard/2.20.2/INSTALL/picard.jar SamToFastq I=${i} FASTQ=fastq_trimmed/${n}_1.fq SECOND_END_FASTQ=fastq_trimmed/${n}_2.fq UNPAIRED_FASTQ=fastq_trimmed/${n}_3.fq; done
