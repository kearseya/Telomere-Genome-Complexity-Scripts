#!/bin/bash

## run dysgu 
ref="/mnt/references/hg38.fa"
cov="python3 dysgu/scripts/coverage2bed.py"

s=(`ls breast/*.bam`)
echo -e "--------------------\n\tLIST (${#s[@]})\n\n${s[*]}\n"
echo -e "--------------------\n\n\tSTART\n\n--------------------"
for i in ${s[*]}
	do
	## run dysgu
	b=$(basename ${i} .bam)
	if [[ -e "dysgu_out/${b}.vcf" ]]; then
		echo "${b} already exists"
		continue
	fi
	echo "Running dysgu ${b}"
	(dysgu run -p4 -v2 --metrics ${ref} tmp_${b} ${i} -o ${b}.vcf 2> ${b}.log ; ${cov} -w tmp_${b} > ${b}_cov.bed ; rm -rf tmp_${b} ; mv ${b}.vcf /mnt/breast/dysgu_out/ ; mv ${b}.log /mnt/breast/dysgu_out/) &
	echo ${b}
	jobs=($(jobs -p))
	while (( ${#jobs[*]} >= 8 ))
	do
		sleep 30
		jobs=($(jobs -p))
	done
done

## separate loop because of basename
s=(`ls breast/*.cram`)
echo -e "--------------------\n\tLIST (${#s[@]})\n\n${s[*]}\n"
echo -e "--------------------\n\n\tSTART\n\n--------------------"
for i in ${s[*]}
	do
	b=$(basename ${i} .cram)
	if [[ -e "dysgu_out/${b}.vcf" ]]; then
		echo "${b} already exists"
		continue
	fi

	echo "Running dysgu ${b}"
	(dysgu run -p4 -v2 --metrics ${ref} tmp_${b} ${i} -o ${b}.vcf 2> ${b}.log ; ${cov} -w tmp_${b} > ${b}_cov.bed ; rm -rf tmp_${b} ; mv ${b}.vcf /mnt/breast/dysgu_out/ ; mv ${b}.log /mnt/breast/dysgu_out/) &
	echo ${b}
	jobs=($(jobs -p))
	while (( ${#jobs[*]} >= 8 ))
	do
		sleep 30
		jobs=($(jobs -p))
	done
done


