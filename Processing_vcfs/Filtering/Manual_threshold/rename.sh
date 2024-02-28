#!/bin/bash

if [[ $# != 3 ]]; then
	echo "Usage: ./rename.sh indir {format}_{string} vcf"
	exit
fi

vars=$(echo $2 | tr "_" "\n")

# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample

function get_index () {
	format_order=$(echo ${line} | awk -F" " '{print $9}')
	IFS=':' read -ra order <<< ${format_order}
	for idx in "${!order[@]}"; do
	   if [[ "${order[${idx}]}" = "${variable}" ]]; then
		   pos="${idx}";
	   fi
	done
}

function return_value () {
	if [[ ${v} == *"{"* ]]; then
		variable=$(echo ${v} | grep -oP "(?<=\{).*?(?=\})")
	else
		val=${v}
	fi
	if [[ "${variable}" = "chrom" ]]; then
		val=$(echo ${line} | awk -F" " '{print $1}')
	elif [[ "${variable}" = "pos" ]]; then
		val=$(echo ${line} | awk -F" " '{print $2}')
	elif [[ "${variable}" = "id" ]]; then
		val=$(echo ${line} | awk -F" " '{print $3}')
	elif [[ "${variable}" = "ref" ]]; then
		val=$(echo ${line} | awk -F" " '{print $4}')
	elif [[ "${variable}" = "alt" ]]; then
		val=$(echo ${line} | awk -F" " '{print $5}')
	elif [[ "${variable}" = "qual" ]]; then
		val=$(echo ${line} | awk -F" " '{print $6}')
	elif [[ "${variable}" = "filter" ]]; then
		val=$(echo ${line} | awk -F" " '{print $7}')
	elif [[ "${variable}" = info* ]]; then
		variable=$(echo ${variable} | grep -oP "(?<=\.).*$")
		val=$(echo ${line} | awk -F" " '{print $8}' | grep -oP "(?<=${variable}\=).*?(?=;)")
	elif [[ "${variable}" = format* ]]; then
		variable=$(echo ${variable} | grep -oP "(?<=\.).*$")
		get_index ${line} ${variable}
		IFS=':' read -ra format_vals <<< $(echo ${line} | awk -F" " '{print $10}')
		val="${format_vals[${pos}]}"
	fi
}

for i in ${1}/*.png; do
	x=$2
	id=$(echo ${i} | grep -oP "(?<=\~)\d+(?=\.png$)")
	line=$(awk -v idv="${id}" -F"\t" '{if ($3 == idv) {print}}' $3)
	for v in ${vars}; do
		val=""
		return_value ${v} ${val} ${line}
		x=$(echo ${x} | sed "s/${v}/${val}/g")
	done
	mv ${i} "icgc/${x}.png"
done
