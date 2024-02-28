#!/bin/bash

export STORAGE_PROFILE=collab
export PATH=$PATH:/home/ubuntu/.local/bin

readarray -t s < /home/ubuntu/normal_ids_test.txt

echo -e "--------------------\n\tLIST (${#s[@]})\n\n${s[*]}\n"
echo -e "--------------------\n\n\tSTART\n\n--------------------"

sc=/home/ubuntu/score-client/bin/score-client

for i in ${s[*]}
	do
	echo "attempting download"
	${sc} download --object-id ${i} --output-dir /home/ubuntu || echo "${i} failed"; continue	
	echo "download successful"
	bam_file=*.bam
	b=$(basename ${bam_file} .bam)
	telseq ${bam_file} -o telseq_out/${b}.telseq
	samtools view -hb --region-file tumours.bed ${bam_file} -o trimmed_normal/${b}.bam
	samtools index trimmed_normal/${b}.bam
	rm *.bam*
done

# ftp transfer trimmed and telseq files after analysis
