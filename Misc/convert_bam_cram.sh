#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH -o /home/User/misc/convert_bam_cram.%J
#SBATCH -e /home/User/misc/convert_bam_cram.%J
#SBATCH --job-name=convert_bam_cram
#SBATCH --account=ProjectID
#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=ALL

module load samtools/1.9

for i in $(ls /scratch/ProjectDir/aligns/*.bam); 
do name=$(basename $i .bam);
samtools view -C -T /scratch/ProjectDir/hg38.fa ${i} > /scratch/ProjectDir/aligns/${name}.cram;
rm ${i};
done
