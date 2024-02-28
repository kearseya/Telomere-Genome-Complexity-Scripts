#!/usr/bin/python3
import multiprocessing
from io import StringIO
from csv import writer
import re
import os
import sys
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
import pysam
import dysgu
from dysgu.python_api import vcf_to_df
import click



def read_from_inputfile(path):
	if path != '-' and isinstance(path, str):
		try:
			with open(path, "r") as fin:
				for line in fin:
					yield line
		except UnicodeDecodeError:
			with gzip.open(path, 'r') as fin:
				for line in fin:
					yield line.decode('ascii')
	elif path == '-':
		for line in sys.stdin:
			yield line
	else:
		for line in path:
			yield line


def combined_vcf_to_df(path):
	#path = "test"
	fin = read_from_inputfile(path)
	# Parse sample
	header = ""
	last = ""
	keep_individual_meta = False
	for line in fin:
		if line[:2] == "##":
			header += line
		else:
			header += "\t".join(line.split("\t")[:9])
			last = line
			break

	samples = last.strip().split("\t")[9:]
	n_samples = len(samples)

	if path == '-':
		path_h = sys.stdin
	else:
		path_h = path

	df = pd.read_csv(path_h, index_col=None, comment="#", sep="\t", header=None)
	df = df.rename({0: "chrA", 1: "posA", 2: "event_id", 3: "ref_seq", 4: "variant_seq", 5: "filter"}, axis=1)

	info = []
	for k in list(df[7]):
		if k:
			info.append(dict(i.split("=") for i in k.split(";") if "=" in i))
		else:
			info.append({i: None for i in info[-1]})

	info_df = pd.DataFrame.from_records(info)

	n_fields = len(df[8][0].split(":"))
	sample_info = []
	for n, s in enumerate(samples):
		for idx, (k1, k2), in enumerate(zip(df[8], df[n+9])):  # Overwrite info column with anything in format
			if n == 0:
				sample_info.append([{i: j for i, j in zip(k1.split(":"), k2.split(":"))}])
			else:
				sample_info[idx].append({i: j for i, j in zip(k1.split(":"), k2.split(":"))})

	sample_df = pd.DataFrame.from_records(sample_info)

	df_cols = ["chrA", "posA", "event_id", "ref_seq", "variant_seq", "filter"]
	info_cols = list(info_df.columns)
	sample_cols = list(df[8][0].split(":"))
	parsed = pd.DataFrame(columns = list(itertools.chain(df_cols, ["sample"], info_cols, sample_cols)))
	n_cols = len(list(itertools.chain(df_cols, ["sample"], info_cols, sample_cols)))
	empty_row = [None]*n_cols
	pidx = 0
	output_csv = StringIO()
	csv_writer = writer(output_csv)
	csv_writer.writerow(list(itertools.chain(df_cols, ["sample"], info_cols, sample_cols)))
	for idx in range(len(df)):
		copies = []
		for i, s in enumerate(samples):
			if sample_info[idx][i]["GT"] != "0":
				copies.append((s, i))
		print(idx, end="\r")
		for s in copies:
			## reliable but slow
			#parsed.loc[pidx] = empty_row
			#parsed.loc[parsed.index == pidx, df_cols] = df.iloc[[idx]][df_cols].set_index(pd.Index([pidx]))
			#parsed.loc[parsed.index == pidx, ["sample"]] = s[0]
			#parsed.loc[parsed.index == pidx, info_cols] = info_df.iloc[[idx]][info_cols].set_index(pd.Index([pidx]))
			#parsed.loc[parsed.index == pidx, sample_cols] = pd.DataFrame.from_records(sample_df[s[1]][idx], index=pd.Index([pidx])) 
			#samples_df.loc[samples_df.index == s][sample_cols].set_index(pd.Index([pidx]))
			## faster but slows down
			#new_row = pd.concat([df.iloc[[idx]][df_cols].set_index(pd.Index([pidx])), 
			#pd.DataFrame({"sample": s[0]}, pd.Index([pidx])),
			#info_df.iloc[[idx]].set_index(pd.Index([pidx])),
			#pd.DataFrame.from_records(sample_df[s[1]][idx], index=pd.Index([pidx]))], axis=1)
			#parsed = pd.concat([parsed, new_row])
			if keep_individual_meta == False:
				new_row = pd.concat([df.iloc[[idx]][df_cols].set_index(pd.Index([pidx])), pd.DataFrame({"sample": "common"}, pd.Index([pidx])), info_df.iloc[[idx]].set_index(pd.Index([pidx]))], axis=1)
				#print(new_row)
			else:
				new_row = pd.concat([df.iloc[[idx]][df_cols].set_index(pd.Index([pidx])), info_df.iloc[[idx]].set_index(pd.Index([pidx])), pd.DataFrame({"sample": s[0]}, pd.Index([pidx])),
				pd.DataFrame.from_records(sample_df[s[1]][idx], index=pd.Index([pidx]))], axis=1)
			if len(new_row.values.flatten().tolist()) >= 20:
				csv_writer.writerow(new_row.values.flatten().tolist())
			else:
				print(new_row)
			del new_row
			#pd.DataFrame({"sample": s[0]}, pd.Index([pidx])),
			#info_df.iloc[[idx]],
			#pd.DataFrame.from_records(sample_df[s[1]][idx], index=pd.Index([pidx]))], axis=1)

			pidx += 1
	#print(output_csv)
	output_csv.seek(0)
	df = pd.read_csv(output_csv)
	#df = parsed.copy()
	#df.to_csv(path+"_pre_drop_conv.csv")
	#df = df[df['chrA'].notna()]

	col_map = {"chrA": ("chrA", str),
			   "posA": ("posA", np.int64),
			   "CHR2": ("chrB", str),
			   "GRP": ("grp_id", np.int64),
			   "NGRP": ("n_in_grp", np.int64),
			   "END": ("posB", np.int64),
			   "CHR2_POS": ("posB_tra", np.int64),
			   "sample": ("sample", str),
			   "event_id": ("event_id", np.int64),
			   "KIND": ("kind", str),
			   "SVTYPE": ("svtype", str),
			   "CT": ("join_type", str),
			   "CIPOS95": ("cipos95A", np.int64),
			   "CIEND95": ("cipos95B", np.int64),
			   "NMP": ("NMpri", float),
			   "NMB": ("NMbase", float),
			   "NMS": ("NMsupp", float),
			   "MAPQP": ("MAPQpri", float),
			   "MAPQS": ("MAPQsupp", float),
			   "NP": ("NP", np.int64),
			   "OL": ("query_overlap", np.int64),
			   "MAS": ("maxASsupp", np.int64),
			   "SU": ("su", np.int64),
			   "WR": ("spanning", np.int64),
			   "PE": ("pe", np.int64),
			   "SR": ("supp", np.int64),
			   "SC": ("sc", np.int64),
			   "BND": ("bnd", np.int64),
			   "SQC": ("sqc", float),
			   "SCW": ("scw", float),
			   "SQR": ("clip_qual_ratio", float),
			   "RT": ("type", str),
			   "BE": ("block_edge", np.int64),
			   "COV": ("raw_reads_10kb", float),
			   "MCOV": ("mcov", float),
			   "LNK": ("linked", np.int64),
			   "CONTIGA": ("contigA", str),
			   "CONTIGB": ("contigB", str),
			   "ref_seq": ("ref_seq", str),
			   "variant_seq": ("variant_seq", str),
			   "GC": ("gc", float),
			   "NEIGH": ("neigh", np.int64),
			   "NEIGH10": ("neigh10kb", np.int64),
			   "REP": ("rep", float),
			   "REPSC": ("rep_sc", float),
			   "LPREC": ("svlen_precise", np.int64),
			   "NEXP": ("n_expansion", np.int64),
			   "STRIDE": ("stride", np.int64),
			   "EXPSEQ": ("exp_seq", str),
			   "RPOLY": ("ref_poly_bases", np.int64),
			   "GT": ("GT", str),
			   "GQ": ("GQ", object),
			   "RB": ("ref_bases", np.int64),
			   "SVLEN": ("svlen", np.int64),
			   "PS": ("plus", np.int64),
			   "MS": ("minus", np.int64),
			   "SBT": ("strand_binom_t", float),
			   "PROB": ("prob", float),
			   "NG": ("n_gaps", float),
			   "NSA": ("n_sa", float),
			   "NXA": ("n_xa", float),
			   "NMU": ("n_unmapped_mates", np.int64),
			   "NDC": ("double_clips", np.int64),
			   "RMS": ("remap_score", np.int64),
			   "RED": ("remap_ed", np.int64),
			   "BCC": ("bad_clip_count", np.int64),
			   "STL": ("n_small_tlen", np.int64),
			   "RAS": ("ras", np.int64),
			   "FAS": ("fas", np.int64),
			   "ICN": ("inner_cn", float),
			   "OCN": ("outer_cn", float),
			   "CMP": ("compress", float),
			   "FCC": ("fcc", float),
			   "RR": ("ref_rep", float),
			   "JIT": ("jitter", float),
			   "LEFT_SVINSSEQ": ("left_ins_seq", str),
			   "RIGHT_SVINSSEQ": ("right_ins_seq", str),
			   }
	df.rename(columns={k: v[0] for k, v in col_map.items()}, inplace=True)
	df["GQ"] = pd.to_numeric(df["GQ"], errors='coerce').fillna(".")
	for k, dtype in col_map.values():
		if k in df:
			if df[k].dtype != dtype:
				if dtype == str:
					df[k] = df[k].fillna("")
				else:
					df[k] = df[k].fillna(0)
				try:
					df[k] = df[k].astype(dtype)
				except ValueError or OverflowError:
					#print(list(df[k][0]))
					raise ValueError("Problem for feature {}, could not intepret as {}".format(k, dtype))
	if "contigA" not in df:
		df["contigA"] = [""] * len(df)
	if "contigB" not in df:
		df["contigB"] = [""] * len(df)

	if 'posB_tra' in df:
		df["posB"] = [i if svt != 'TRA' else j for i, j, svt in zip(df['posB'], df['posB_tra'], df['svtype'])]
		del df['posB_tra']
	df["posA"] = df["posA"].astype(int) - 1  # convert to 0-indexing
	df["posB"] = df["posB"].astype(int) - 1

	write = True
	if write == True:
		df.to_csv(path+"_conv.csv")

	#print(df)
	return df, header, n_fields



@click.command()
@click.argument("filename", nargs=1)
@click.argument("merged", nargs=1)
@click.option("-o", "--outdir", default=".",
		help="Output direectory (default current)")
@click.option("-f", "--fmt", type=click.Choice(["csv", "vcf"]), default = "vcf",
		help="Output format (vcf or csv)")
@click.option("-p", "--plot", is_flag=True, default=False, show_default=True,
		help="Remove supplementary matertial")
@click.option("-s", "--suffix", default="_filtered",
		help="Output file suffix")
@click.option("-c", "--count", default=1,
		help="Threshold for number of sample to keep variant")
@click.option("-m", "--low-mem", is_flag=True, default=False, show_default=True,
		help="Don't store all variants from one parse, rather parse for each chromosome")
@click.option("--region-based", is_flag=True, default=False, show_default=True,
		help="Only use proximity of mutations to filter")
def remove_common_variants(filename, merged, outdir, fmt, plot, suffix, count, low_mem, region_based):
	"""Removes variants from vcf found in a merged vcf (common variants)"""
	if region_based == False:
		sample = str(os.path.basename(filename).rsplit(".", 1)[0]) 
		out_file_path = os.path.join(outdir, sample+suffix+"."+fmt)

		if count > 1:
			print("Reading file...")
			filtered = pd.read_csv(merged, comment="#", sep="\t", index_col=None, header=None)
			counts = []
			cols = filtered.columns[9:-1]
			print("Counting variant frequency...")
			total_vars = len(filtered)
			for index, row in filtered.iterrows():
				rc = 0
				for col in cols:
					if re.search('^0[:0]*', row[col]) == None:
						rc += 1
				counts.append(rc)
				print(f"{index}/{total_vars}", end="\r")
			del filtered
			
			if plot == True:
				print("Plotting Variants...")
				fig, (ax1, ax2) = plt.subplots(2)
				ax1.hist(counts)
				ax1.axvline(count, color='k', linestyle='dashed', linewidth=1)
				ax2.hist([i for i in counts if i >= count])
				plt.show()
				plotted_count = input("Revised count: ")
				if plotted_count == "":
					plotted_count = count
				if int(count) != int(plotted_count):
					count = int(plotted_count)

			print("Getting header index...")
			filtered_file_base = str(os.path.basename(merged).rsplit(".", 1)[0])
			filtered_file = os.path.join(outdir, filtered_file_base+f"_filtered_{count}.vcf")
			with open(merged, "r") as f:
				for line in f:
					if line.startswith("#"):
						counts.insert(0, count+1)
					else:
						break
			
			print("Filtering file for low frequency variants...")
			filtered_vars = 0
			with open(merged, "r") as f:
				with open(filtered_file, "w") as o:
					for idx, line in enumerate(f):
						if counts[idx] >= count:
							o.write(line)
							filtered_vars += 1
						print(f"{filtered_vars}/{total_vars}", end="\r")
			print("")
			print(f"Removed {total_vars-filtered_vars} variants, writted to {filtered_file}")

			merged = filtered_file


		sample_fmt = filename[-3:]
		if sample_fmt == "csv":
			in_df = pd.DataFrame.from_csv(filename)
		elif sample_fmt == "vcf":
			in_df, in_header, in_fields = vcf_to_df(filename)
		else:
			raise ValueError(f"{sample_fmt} invalid file type")

		merged_fmt = merged[-3:]
		if merged_fmt == "csv":
			merged_df = pd.read_csv(merged)
		elif merged_fmt == "vcf":
			merged_df, merged_header, merged_fields = combined_vcf_to_df(merged)
		else:
			raise ValueError(f"{merged_fmt} invalid file type")

		#print(in_df)
		#print(merged_df)
		all_df = dysgu.merge_dysgu_df(in_df, merged_df)
		print(all_df)
		all_df.to_csv("all.csv")
		

	if region_based == True:
		## create new file with sample header
		sample = str(os.path.basename(filename).rsplit(".", 1)[0]) 
		out_file_path = os.path.join(outdir, sample+suffix+".vcf")
		#if os.path.exists(out_file_path):
		#		print("Outfile exists, will not overwrite")
		#		quit()
		in_file = pysam.VariantFile(filename, "r")
		out_file = pysam.VariantFile(out_file_path, "w", header=in_file.header)
		#infile.close()
		merged_file = pysam.VariantFile(merged, "r")

		## create dictionaries for locations for common and tumour variant locations
		common, tumour = {}, {}
		for contig in list(in_file.header.contigs):
			common[contig], tumour[contig] = [], []
		if low_mem == False:	
			## iterate tumour file for variant positions
			total = 0
			for variant in in_file:
					tumour[variant.chrom].append((variant.start, variant.stop, variant.id))
					total += 1
			
			## iterate merged file for common variants
			total_common = 0
			for variant in merged_file:
				common[variant.chrom].append((variant.start, variant.stop))
				total_common += 1
			merged_file.close()

		
		## check overlap with ids
		unique_ids = set()
		for contig in tumour.keys():
			if low_mem == True:
				common, tumour = {}, {}
				common["contig"], tumour["contig"] = [], []
				for variant in merged_file.fetch(contig=contig):
					common["contig"].append((variant.start, variant.stop, -1))
				for variant in in_file.fetch(contig=contig):
					tumour["contig"].append((variant.start, variant.stop, variant.id))
				contig = "contig"
				
			intervals = []
			for region in tumour[contig]:
				intervals.append((contig, region[0], region[1], region[2]))
			for region in common[contig]:
				intervals.append((contig, region[0], region[1], -1))
			overlap = dysgu.merge_intervals(intervals, add_indexes=True)
			if len(overlap) > 0:
				for i in overlap:
					if -1 not in i[3]:
						unique_ids.add(i[3][0])
		print("size of tumour: ", sys.getsizeof(tumour))
		print("size of normal: ", sys.getsizeof(common))
		print(f"{len(unique_ids)}/{total} reads from tumour unique from {total_common} common variants")		
		## write unique reads to file
		in_file.reset()
		for variant in in_file:
			if variant.id in unique_ids:
				out_file.write(variant)
				unique_ids.remove(variant.id)
		
		in_file.close()
		out_file.close()


if __name__ == "__main__":
	remove_common_variants()
	#procs = 1	# Number of processes to create

	## Create a list of jobs and then iterate through
	## the number of processes appending each process to
	## the job list 
	#jobs = []
	#for i in range(0, procs):
	#	out_list = list()
	#	process = multiprocessing.Process(target=remove_common_variants, 
	#							  args=())
	#	jobs.append(process)

	## Start the processes (i.e. calculate the random number lists)		
	#for j in jobs:
	#	j.start()

	## Ensure all of the processes have finished
	#for j in jobs:
	#	j.join()
