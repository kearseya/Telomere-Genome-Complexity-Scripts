#!/bin/bash

module load jellyfish

jellyfish count -m 21 -s 3G -L 1 -U 1 /scratch/ProjectDir/reference_genomes/hg38.fa -o hg38_21.jf
jellyfish dump hg38_21.jf -o filter.fa
