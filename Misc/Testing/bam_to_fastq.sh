#!/bin/bash

module load picard/2.20.2

for i in $(ls aligns/*bam | head); do n=$(basename ${i} .bam); java -jar /software/genomics/picard/2.20.2/INSTALL/picard.jar SamToFastq I=${i} FASTQ=fastq_db/${n}_1.fq SECOND_END_FASTQ=fastq_db/${n}_2.fq UNPAIRED_FASTQ=fastq_db/${n}_3.fq; done

