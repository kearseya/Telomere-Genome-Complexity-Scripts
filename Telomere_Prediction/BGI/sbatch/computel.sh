#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH -o /home/User/telomere_prediction/telomerecat_std_logs/out.%J
#SBATCH -e /home/User/telomere_prediction/telomerecat_std_logs/err.%J
#SBATCH --job-name=telomerecat
#SBATCH --account=ProjectID
#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=ALL

computel_path=/scratch/User/miniconda3/bin
data_path=/scratch/ProjectDir/hap1/
out_path=/scratch/ProjectDir/output/telomerecat/

for i in $(basename -s .bam ${data_path}*.bam | uniq); 
do 
	${computel_path}/telomerecat bam2length -p 10 --output ${out_path} ${i}.bam --temp_dir /scratch/ProjectDir/tmp 
done
