import os
import pysam

in_dir = "/scratch/ProjectDir/hap1/"

bam_files = os.listdir(in_dir)
bam_files = [f for f in bam_files if f.endswith(".bam")]

def index_stats(file_in, sample, in_dir):
    f = pysam.AlignmentFile(in_dir+file_in)
    s = f.get_index_statistics()
    total = 0
    ref_len = 0
    for i in s:
        total += i.mapped
        ref_len += f.get_reference_length(i.contig)
    cov = total*98/ref_len
    print(sample,cov)
    return cov

for file in bam_files:
    fn = str.split(file, "/")[-1]
    sample = str.split(fn, ".")[0]
    index_stats(file, sample, in_dir)
