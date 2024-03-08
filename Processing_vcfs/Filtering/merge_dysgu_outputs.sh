#! /bin/bash

### before dysgu filter was made

merge_individual_tumour_normal=false
if [ "$merge_individual_tumour_normal" = true ]; then
while IFS="," read -r normal tumour
do
	echo "MERGING normal: ${normal}, tumour ${tumour}"
	dysgu merge -o ${normal}_${tumour}.csv -f csv breast/dysgu_out/${normal}.vcf breast/dysgu_out/${tumour}.vcf
	mv ${normal}_${tumour}.csv merged_csv
done < <(cut -d "," -f1,2 sample_pairs.csv | tail -n +2)
fi

# merge tumour and normal
merge_all_samples=false
if [ "$merge_all_samples" = true ]; then
echo "MERGING ALL"
dysgu merge -o all_merged.vcf breast/dysgu_out/*.vcf
fi

# drop common variants
merge_against_all=false
if [ "${merge_against_all}" = true ]; then
for i in breast/dysgu_out/*.vcf;
do
	b=$(basename ${i} .vcf)
	echo "ALL MERGE: ${b}"
	awk 'FNR==NR {a[$0]++; next} !($0 in a)' all_merged.vcf ${i} > all_merge_drop/${b}_drop_common.vcf
done
fi

# merge all tumour samples
merge_all_tumours=false
if [ "$merge_all_tumours" = true ]; then
tumour_files=()
while IFS="," read -r normal tumour
do
	if [ -f "all_merge_drop/${tumour}_drop_common.vcf" ]; then 
		tumour_files+=( "all_merge_drop/${tumour}_drop_common.vcf" )
	fi
done < <(cut -d "," -f1,2 sample_pairs.csv | tail -n +2)

echo "MERGING TUMOURS"
dysgu merge -o tumour_merged.vcf ${tumour_files[*]}
fi

# drop common variants
merge_against_tumours=true
if [ "$merge_against_tumours" = true ]; then
tumour_files=()
while IFS="," read -r normal tumour
do
	if [ -f "all_merge_drop/${tumour}_drop_common.vcf" ]; then 
		tumour_files+=( "all_merge_drop/${tumour}_drop_common.vcf" )
	fi
done < <(cut -d "," -f1,2 sample_pairs.csv | tail -n +2)

for i in ${tumour_files[*]}
do
	b=$(basename ${i} _drop_common.vcf)
	echo "TUMOUR MERGE: ${b}"
	awk 'FNR==NR {a[$0]++; next} !($0 in a)' tumour_merged.vcf ${i} > tumour_merge_drop/${b}_drop_common_2.vcf
done
fi


# run script
run_filter=false
if [ "$run_filter" = true ]; then
python3 filter_reads_for_svs.py
fi



