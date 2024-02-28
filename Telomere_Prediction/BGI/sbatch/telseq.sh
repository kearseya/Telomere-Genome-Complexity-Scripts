#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH -o /home/User/telomere_prediction/predictors/time_test/telseq_out.%J
#SBATCH -e /home/User/telomere_prediction/predictors/time_test/telseq_err.%J
#SBATCH --job-name=telseq
#SBATCH --account=ProjectID
#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=END

telseq_path=/scratch/User/miniconda3/bin
data_path=/scratch/ProjectDir/aligns/
out_path=/scratch/ProjectDir/output/telseq

${telseq_path}/telseq -f testlist.txt > ${out_path}/telseq_hap1.csv

