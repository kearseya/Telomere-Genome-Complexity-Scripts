#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH -o /home/User/telomere_prediction/predictors/time_test/telomerecat_out.%J
#SBATCH -e /home/User/telomere_prediction/predictors/time_test/telomerecat_err.%J
#SBATCH --job-name=telomerecat
# #SBATCH --account=ProjectID
#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=END

telcat_path=/scratch/User/miniconda3/bin
data_path=/scratch/ProjectDir/aligns/
out_path=/scratch/ProjectDir/output/telomerecat/

cd ${data_path}

# for i in $(basename -s .bam ${data_path}*.bam | uniq); 
# do 
# ${telcat_path}/telomerecat bam2length -p 10 --output ${out_path} ${i}.bam --temp_dir /scratch/ProjectDir/tmp 
# done
# 

i=DB208

${telcat_path}/telomerecat bam2length -p 10 --output ${out_path} ${i}.bam --temp_dir /scratch/ProjectDir/tmp 

