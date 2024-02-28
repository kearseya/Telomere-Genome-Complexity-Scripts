#!/bin/bash
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8

module load jellyfish/2.2.10
module load samtools
module load python
module load htslib

for i in /scratch/ProjectDir/aligns/*.bam; 
    do b=$(basename ${i} .bam); 
    if [ ! -f /scratch/ProjectDir/jellyfish/done/${b}.bin.done ];
    then
        jellyfish count -m 21 -s 3G -t 8 --if filter.fa -o /scratch/ProjectDir/jellyfish/${b}.jf <(samtools fasta ${i}); 
        jellyfish dump -c ${b}.jf | python convert.py ${b}.bin;
    fi
    rm ${b}.jf;
done

# jellyfish count -m 13 -s 1G -t 8 --sam ${i} -o /scratch/ProjectDir/jellyfish/${b}.jf; 

for i in /scratch/ProjectDir/aligns/*.cram; 
    do b=$(basename ${i} .cram); 
    if [ ! -f /scratch/ProjectDir/jellyfish/done/${b}.bin.done ];
    then
        jellyfish count -m 21 -s 3G -t 8 --if filter.fa -o /scratch/ProjectDir/jellyfish/${b}.jf <(samtools fasta ${i}); 
        jellyfish dump -c ${b}.jf | python convert.py ${b}.bin;
        rm ${b}.jf;
    fi
done

