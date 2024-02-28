import numpy as np
from numpy.random import choice
import pysam
import os
import pandas as pd

in_dir = "/scratch/ProjectDir/aligns/"
out_dir = "/scratch/ProjectDir/output/pysam_coverage/"
cov_regions = pd.read_csv("hg38_GC_0.48_0.52_bin_1000_1000").to_dict()

coverages = {}
chr_list = ["chr"+str(i) for i in range(1, 22)]+["chrX", "chrY"]
num_chr = len(chr_list)

bam_files = os.listdir(in_dir)
bam_files = [f for f in bam_files if f.endswith(".bam")]
total_files = len(bam_files)

print(bam_files)


## Random
coverages = {}
chr_list = ["chr"+str(i) for i in range(1, 22)]+["chrX", "chrY"]
num_chr = len(chr_list)

if os.path.isdir(in_dir) == True:
    bam_files = os.listdir(in_dir)
    bam_files = [f for f in bam_files if f.endswith(".bam")]
    total_files = len(bam_files)
    for num, bam_file in enumerate(bam_files):
        print("Counting coverage for:   ", num, "/", total_files, "     (", bam_file, ")            ", end="\r")
        bamfile = pysam.AlignmentFile(os.path.join(in_dir, bam_file), "rb")
        basename = os.path.basename(bam_file).split(".")[0]
        header = bamfile.header
        head = {}
        for c in header.references:
            head[c] = header.get_reference_length(c)

        chr_lens = [int(head[c]) for c in chr_list]
        total_lens = sum(chr_lens)
        chr_lens = [l/total_lens for l in chr_lens]
        #print(chr_lens)
        #print(sum(chr_lens))
        coverages[basename] = []
        for x in range(10000):
            chrom = choice(chr_list, p=chr_lens)
            limit = head[chrom]
            rand_pos_1 = np.random.randint(0, limit)
            if limit-rand_pos_1 < 1000:
                rand_pos_1 = limit-1000
                rand_pos_2 = limit
            else:
                rand_pos_2 = rand_pos_1 - 1000
            cov = np.sum(bamfile.count_coverage(chrom, rand_pos_1, rand_pos_2))
            coverages[basename].append(cov)
#print(coverages)
coverages = pd.DataFrame(coverages)
coverages.to_csv(os.path.join(out_dir, "10k_aligns_coverages.csv"), index=False)












#for num, bam_file in enumerate(bam_files):
#    print("Counting coverage for:   ", num, "/", total_files, "     (", bam_file, ")", end="\r")
#    bamfile = pysam.AlignmentFile(os.path.join(in_dir, bam_file), "rb")
#    basename = os.path.basename(bam_file).split(".")[0]
#    coverages[basename] = []
#    for x in range(len(cov_regions["chromo"])):
#        coverages[basename].append(np.sum(bamfile.count_coverage(cov_regions["chromo"][x], cov_regions["lower"][x], cov_regions["upper"][x])))
#coverages = pd.DataFrame(coverages)
#coverages.to_csv(os.path.join(out_dir, "hap1_gc_coverages.csv"), index=False)
#



#for num, bam_file in enumerate(bam_files):
#    print("Counting coverage for:   ", num, "/", total_files, "     (", bam_file, ")", end="\r")
#    bamfile = pysam.AlignmentFile(os.path.join(in_dir, bam_file), "rb")
#    basename = os.path.basename(bam_file).split(".")[0]
#    header = pysam.view("-H", os.path.join(in_dir, bam_file)).split("\n")
#    head = {}
#    for line in header:
#        if line.startswith("@SQ") == True:
#            cut = line.split("SN:")[1]
#            chrom = cut.split("LN:")[0].strip()
#            ln = cut.split("LN:")[1].strip()
#            head[chrom] = ln
#
#    chr_lens = [int(head[c]) for c in chr_list]
#    total_lens = sum(chr_lens)
#    chr_lens = [l/total_lens for l in chr_lens]
#    coverages[basename] = []
#    for x in range(1000):
#        chrom = choice(chr_list, p=chr_lens)
#        limit = int(head[chrom])
#        rand_pos_1 = np.random.randint(0, limit)
#        if limit-rand_pos_1 < 1000:
#            rand_pos_1 = limit-1000
#            rand_pos_2 = limit
#        else:
#            rand_pos_2 = rand_pos_1 + 1000
#        cov = np.sum(bamfile.count_coverage(chrom, rand_pos_1, rand_pos_2))
#        coverages[basename].append(cov)
#    print(sum(coverages[basename])/1000)
#print(coverages)
#coverages = pd.DataFrame(coverages)
#coverages.to_csv(os.path.join(out_dir, "coverages.csv"), index=False)
