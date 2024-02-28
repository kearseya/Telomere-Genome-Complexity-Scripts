#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH -o /home/User/telomere_prediction/predictors/time_test/qmotif_out.%J
#SBATCH -e /home/User/telomere_prediction/predictors/time_test/qmotif_err.%J
#SBATCH --job-name=qmotif
#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=END

qmotif_path=/home/User/tools/adamajava/qmotif/build/flat/
data_path=/scratch/ProjectDir/aligns/
out_path=/scratch/ProjectDir/output/qmotif/time_test

for i in $(basename -s .bam ${data_path}*.bam | uniq); 
do 
java -Xmx20g -jar ${qmotif_path}qmotif.jar\
	-n 8 \
	-bam ${data_path}${i}.bam \
	-bai ${data_path}${i}.bam.bai \
	--log ${out_path}${i}.qmotif.log \
	-ini /home/User/tools/adamajava/qmotif/qmotif.ini \
	-o ${out_path}${i}.qmotif.xml \
	-o ${out_path}${i}.telomere.bam; 
done

# i=DB208

# java -Xmx20g -jar ${qmotif_path}qmotif.jar\
# 	-n 8 \
# 	-bam ${data_path}${i}.bam \
# 	-bai ${data_path}${i}.bam.bai \
# 	--log ${out_path}${i}.qmotif.log \
# 	-ini /home/User/tools/adamajava/qmotif/qmotif.ini \
# 	-o ${out_path}${i}.qmotif.xml \
# 	-o ${out_path}${i}.telomere.bam;
