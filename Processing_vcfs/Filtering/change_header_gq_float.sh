#!/bin/bash

for i in ${1}/*.vcf;
do
	sed -i 's/GQ,Number=1,Type=Integer/GQ,Number=1,Type=Float/g' ${i}
done
