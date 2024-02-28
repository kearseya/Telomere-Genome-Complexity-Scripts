#!/bin/bash

for i in ${1}/*.vcf;
do
	sed -e "88r header_contigs.txt" ${i} > ${2}/$(basename ${i})
done
