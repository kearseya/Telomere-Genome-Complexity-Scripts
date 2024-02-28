#!/bin/bash

module load bwa

bwa fastmap ../reference_genomes/hg19.fa filter.fa | cut -s -d$'\t' -f5 | sed '/^$/d' > hg19_mappings

