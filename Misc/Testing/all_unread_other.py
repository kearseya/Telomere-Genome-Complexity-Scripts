import pysam
import pandas as pd
import os
from collections import deque

in_dir = "/scratch/ProjectDir/aligns/"
coords = pd.read_csv("hg38_all_scan_telomeres.tsv", sep="\t")
out_dir = "/scratch/ProjectDir/output/pysam_regions/all_unmapped"



def count_variant_repeats(seq, targets, targets_to_array, direction, k=6):
    seq = seq.upper()
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    kmers_c = 0
    blocks = []
    a = [0] * (len(targets_to_array) + 1)  # last column is for 'others'
    for item in kmers:
        if item in targets:
            a[targets[item]] += 1
            kmers_c += 1
        else:
            a[-1] += 1
    i = 0
    while i < len(kmers):
        current = i
        mm = 0
        while current < len(kmers):
            if kmers[current] in targets:
                current += 1
            else:
                break
        if i != current:
            blocks.append(current - i + k - 1)
            i = current + 1
        i += 1

    kmers_7 = [seq[i:i+7] for i in range(len(seq)-7)]
    if direction == "forward":
        counts_7 = len([1 for i in kmers_7 if i == "CCCTTAA"])
        kmers_c += counts_7
        a[targets["CCCTTAA"]] = counts_7

    if direction == "reverse":
        counts_7 = len([1 for i in kmers_7 if i == "TTAAGGG"])
        kmers_c += counts_7
        a[targets["TTAAGGG"]] = counts_7

    if direction == "both":
        counts_7 = len([1 for i in kmers_7 if i == "CCCTTAA"])
        kmers_c += counts_7
        a[targets["CCCTTAA"]] = counts_7
        counts_7 = len([1 for i in kmers_7 if i == "TTAAGGG"])
        kmers_c += counts_7
        a[targets["TTAAGGG"]] = counts_7

    return kmers_c, a, blocks

def tel_tokens(telomere_f, key_index):
    d = deque(telomere_f)
    f_rotations = {"".join(d): key_index}
    for i in range(len(telomere_f)-1):
        d.rotate()
        f_rotations["".join(d)] = key_index
    return f_rotations

def make_rotation_keys(tta):
    variant_rotations = {}
    for key, key_index in tta.items():
        variant_rotations.update(tel_tokens(key, key_index))
    return variant_rotations



def add(a, arr, low, mapq):
    if a.flag & 1284 or a.cigartuples is None: # or a.mapq == 0:
        return arr
    if int(a.reference_start)-low+100 > len(arr)-1:
        return arr
    if mapq == False:
        index_start = a.reference_start-low+100
    if mapq == True:
        index_start = 100
    index_bin = index_start
    for opp, length in a.cigartuples:
        if index_bin > len(arr) - 1:
            break
        if opp == 2:
            index_start += length
            index_bin = index_start
        elif opp == 0 or opp == 7 or opp == 8: # cigar is mapped
            arr[index_bin] += 1
            index_start += length
            index_bin = index_start
            if index_bin < len(arr):
                arr[index_bin] -= 1
    return arr




targets_to_array_f = {"CCCTAA": 0, "CCCTGA": 1, "CCCGAA": 2, "CCCTAC": 3, "CCCTCA": 4, "CCCCAA": 5, "CCCTTA": 6,
                   "CCCTAT": 7, "CCCTAG": 8, "CCCAAA": 9}
targets_to_array_f.update({"CCCTTAA": 10, "CCCACT": 11, "CCCCAT": 12, "CCCGCA": 13, "CCCGCT": 14, "CCCTCT": 15})

targets_to_array_r = {"TTAGGG": 0, "TCAGGG": 1, "TTCGGG": 2, "GTAGGG": 3, "TGAGGG": 4, "TTGGGG": 5, "TAAGGG": 6,
                   "ATAGGG": 7, "CTAGGG": 8, "TTTGGG": 9}
targets_to_array_r.update({"TTAAGGG": 10, "AGTGGG": 11, "ATGGGG": 12, "TGCGGG": 13, "AGCGGG": 14, "AGAGGG": 15})

targets_to_array_both = {}
targets_to_array_both.update(targets_to_array_f)
targets_to_array_both.update(targets_to_array_r)

targets_f = make_rotation_keys(targets_to_array_f)
targets_r = make_rotation_keys(targets_to_array_r)
targets_both = make_rotation_keys(targets_to_array_both)

targets_dict = {"forward": targets_f, "reverse": targets_r}
targets_to_array_dict = {"forward": targets_to_array_f, "reverse": targets_to_array_r}



bam_files = os.listdir(in_dir)
bam_files = [f for f in bam_files if f.endswith(".bam")]
total_files = len(bam_files)
for num, bam_file in enumerate(bam_files):
    print("Trimming:    ", num, "/", total_files, "     (", bam_file,")         ", end="\r")
    bamfile = pysam.AlignmentFile(os.path.join(in_dir, bam_file), "rb")
    tmpfile = pysam.AlignmentFile(os.path.join(out_dir, os.path.basename(bam_file).split(".")[0]+"_tel.bam"), "wb", template=bamfile)
    #chr_list = ["chr"+str(i) for i in range(1, 22)]+["chrX", "chrY"]
    chr_list = list(set(sorted(coords.chrom)))
    for chromo in chr_list:
        chrom_reg = coords[coords["chrom"] == chromo]
        #print(chrom_reg)
        for i in range(len(chrom_reg["chrom"])):
            low = int(chrom_reg["chromStart"].iloc[i])
            up = int(chrom_reg["chromEnd"].iloc[i])
            #print("up: ", up, "     low: ", low)
            for read in bamfile.fetch(chromo, low, up):
                tmpfile.write(read)
    # unmapped reads
    bamfile = pysam.AlignmentFile(os.path.join(in_dir, bam_file), "rb")
    for read in bamfile:#bamfile.fetch("chr1", 0, until_eof=True):
        if read.is_unmapped == True:
            cf, af, bf, = count_variant_repeats(read.query_alignment_sequence.upper(), targets_f, targets_to_array_f, "forward")
            cr, ar, br, = count_variant_repeats(read.query_alignment_sequence.upper(), targets_r, targets_to_array_r, "reverse")
            if cf >= cr:
                if sum(bf)/read.query_alignment_length > 0.5:
                    tmpfile.write(read)
            if cr > cf:
                if sum(br)/read.query_alignment_length > 0.5:
                    tmpfile.write(read)
    bamfile.close()
    tmpfile.close()
