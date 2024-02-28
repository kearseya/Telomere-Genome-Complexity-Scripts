#!/bin/bash

for i in {1..43}; do
	telomerehunter -ibt telomerehunter_workdir/${i}/tumor_TelomerCnt_${i}/${i}_filtered.bam -ibc telomerehunter_workdir/${i}/control_TelomerCnt_${i}/${i}_filtered.bam -o telomerehunter_filtered -p ${i}
done
