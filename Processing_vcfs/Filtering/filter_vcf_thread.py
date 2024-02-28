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
#from dysgu import DysguSV
from dysgu.python_api import vcf_to_df
from dysgu.io_funcs import col_names
import click
import logging
import pkg_resources
from collections import defaultdict
import sortedcontainers
import subprocess
import time
from datetime import timedelta


def check_header_len(fname):
	lines = 0
	with open(fname, "r") as f:
		for line in f:
			if line.startswith("#"):
				lines += 1
			else:
				return lines

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


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
				new_row = pd.concat([df.iloc[[idx]][df_cols].set_index(pd.Index([pidx])), pd.DataFrame({"sample": "common"}, pd.Index([pidx])), info_df.iloc[[idx]].set_index(pd.Index([pidx])), 
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
					print(list(df[k]))
					raise ValueError("Problem for feature {}, could not intepret as {}".format(k, dtype))
	if "contigA" not in df:
		df["contigA"] = [""] * len(df)
	if "contigB" not in df:
		df["contigB"] = [""] * len(df)

	if 'posB_tra' in df:
		df["posB"] = [i if svt != 'TRA' else j for i, j, svt in zip(df['posB'], df['posB_tra'], df['svtype'])]
		del df['posB_tra']
	df["posA"] = df["posA"].astype(int) - 1 # convert to 0-indexing
	df["posB"] = df["posB"].astype(int) - 1

	write = False
	if write == True:
		df.to_csv(path+"_conv.csv")

	#print(df)
	return df, header, n_fields, idx

def make_main_record(r, version, index, format_f, df_rows, add_kind, extended, small_output):
	rep, repsc, lenprec = 0, 0, 1
	mean_prob, max_prob = None, None
	#print("make df_rows", df_rows)
	if len(format_f) > 1:
		#print(format_f)
		#print(df_rows)
		#print(sorted([(int(v["su"]), k) for k, v in df_rows.items()], reverse=True))
		best = sorted([(int(v["su"]), k) for k, v in df_rows.items()], reverse=True)[0][1]
		probs = [v["prob"] for k, v in df_rows.items()]
		mean_prob = np.mean(probs)
		max_prob = np.max(probs)
		r = df_rows[best]
		gc = round(r["gc"], 2)
		if gc == -1:
			gc = "."
		if not small_output:
			rep = r["rep"]
			repsc = r["rep_sc"]
		lenprec = r["svlen_precise"]
		n_expansion = r["n_expansion"]
		stride = r["stride"]
		exp_seq = r["exp_seq"]
		ref_poly = r["ref_poly_bases"]
		overlaps = r["query_overlap"]

		su, pe, sr, sc, bnd, wr = 0, 0, 0, 0, 0, 0
		for row in df_rows.values():
			pe += row["pe"]
			sr += row["supp"]
			sc += row["sc"]
			su += row["su"]
			bnd += row["bnd"]
			wr += row["spanning"]

	else:
		pe = r["pe"]
		sr = r["supp"]
		sc = r["sc"]
		su = r["su"]
		bnd = r["bnd"]
		wr = r["spanning"]
		gc = round(r["gc"], 2)
		if gc == -1:
			gc = "."
		if not small_output:
			rep = r["rep"]
			repsc = r["rep_sc"]
		lenprec = r["svlen_precise"]
		n_expansion = r["n_expansion"]
		stride = r["stride"]
		exp_seq = r["exp_seq"]
		ref_poly = r["ref_poly_bases"]
		overlaps = r["query_overlap"]

	samp = r["sample"]
	read_kind = r["type"]

	r["posA"] = max(1, r["posA"] + 1)  # convert to 1 based indexing
	r["posB"] = max(1, r["posB"] + 1)

	info_extras = []
	if r["chrA"] == r["chrB"]:
		if 'svlen' in r:
			info_extras.append(f"SVLEN={r['svlen']}")
		else:
			info_extras.append(f"SVLEN=0")
	#print(r["contigA"])
	if r["contigA"]:
		info_extras.append(f"CONTIGA={r['contigA']}")
	if r["contigB"]:
		info_extras.append(f"CONTIGB={r['contigB']}")

	if not r["variant_seq"] or r["variant_seq"][0] == "<":
		if "left_ins_seq" in r and r["left_ins_seq"]:
			info_extras.append(f"LEFT_SVINSSEQ={r['left_ins_seq']}")
		if "right_ins_seq" in r and r["right_ins_seq"]:
			info_extras.append(f"RIGHT_SVINSSEQ={r['right_ins_seq']}")

	if add_kind:
		info_extras += [f"KIND={r['kind']}"]

	if not small_output:

		try:
			info_extras.append(f"REP={'%.3f' % float(rep)}")
		except ValueError:	# rep is "."
			info_extras.append("REP=0")

		try:
			info_extras.append(f"REPSC={'%.3f' % float(repsc)}")
		except ValueError:
			info_extras.append("REPSC=0")

	info_extras += [
					f"GC={gc}",
					f"NEXP={n_expansion}",
					f"STRIDE={stride}",
					f"EXPSEQ={exp_seq}",
					f"RPOLY={ref_poly}",
					f"OL={overlaps}",
					f"SU={su}",
					f"WR={wr}",
					f"PE={pe}",
					f"SR={sr}",
					f"SC={sc}",
					f"BND={bnd}",
					f"LPREC={lenprec}",
					f"RT={read_kind}"]

	if r['chrA'] != r['chrB']:
		chr2_pos = f";CHR2_POS={r['posB']}"
	else:
		chr2_pos = ""

	if mean_prob is not None:
		info_extras += [f"MeanPROB={round(mean_prob, 3)}", f"MaxPROB={round(max_prob, 3)}"]

	if small_output:
		fmt_keys = "GT:GQ:MAPQP:SU:WR:PE:SR:SC:BND:COV:NEIGH10:PS:MS:RMS:RED:BCC:FCC:ICN:OCN:PROB"
	elif extended:
		fmt_keys = "GT:GQ:DP:DN:DAP:DAS:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BND:SQC:SCW:SQR:BE:COV:MCOV:LNK:NEIGH:NEIGH10:RB:PS:MS:SBT:NG:NSA:NXA:NMU:NDC:RMS:RED:BCC:FCC:STL:RAS:FAS:ICN:OCN:CMP:RR:JIT:PROB"
	else:
		fmt_keys = "GT:GQ:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BND:SQC:SCW:SQR:BE:COV:MCOV:LNK:NEIGH:NEIGH10:RB:PS:MS:SBT:NG:NSA:NXA:NMU:NDC:RMS:RED:BCC:FCC:STL:RAS:FAS:ICN:OCN:CMP:RR:JIT:PROB"

	if "variant_seq" in r and isinstance(r["variant_seq"], str):
		if r['svtype'] == "INS":
			alt_field = r.variant_seq.upper()

		else:
			alt_field = f"<{r['svtype']}>"
	else:
		alt_field = f"<{r['svtype']}>"

	if "ref_seq" in r and isinstance(r["ref_seq"], str):
		ref_field = r["ref_seq"]
	else:
		ref_field = "."

	rec = [r["chrA"], r["posA"], r["event_id"],
		   ref_field,
		   alt_field,
		   ".", "." if "filter" not in r else r['filter'],
		   # INFO line
		   ";".join([f"SVMETHOD=DYSGUv{version}",
				   f"SVTYPE={r['svtype']}",
				   f"END={r['posB']}" if r['chrA'] == r['chrB'] else f"END={r['posA'] + 1}",
				   f"CHR2={r['chrB']}" + chr2_pos,
				   f"GRP={r['grp_id']}",
				   f"NGRP={r['n_in_grp']}",
				   f"CT={r['join_type']}",
				   f"CIPOS95={r['cipos95A']}",
				   f"CIEND95={r['cipos95B']}",

				   ] + info_extras),
		   fmt_keys  # :PROB
		   ]
	#print("make rec", rec)
	# FORMAT line(s)
	for item in format_f.values():
		rec.append(":".join(map(str, item)))
	
	return rec



def get_fmt(r, extended=True, small_output=False):
	if small_output:
		v = [r["GT"], r["GQ"], r['MAPQpri'], r['su'], r['spanning'], r['pe'], r['supp'], r['sc'], r['bnd'],
			 r['raw_reads_10kb'], r['neigh10kb'], r["plus"], r["minus"], r["remap_score"], r["remap_ed"],
			 r["bad_clip_count"], round(r["fcc"], 3), round(r["inner_cn"], 3), round(r["outer_cn"], 3), r['prob']
			 ]
		return v

	elif extended:
		v = [r["GT"], r['GQ'], r['DP'], r['DN'], r['DApri'], r['DAsupp'], r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
			 r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
			 r['sc'], r['bnd'], round(r['sqc'], 2), round(r['scw'], 1), round(r['clip_qual_ratio'], 3), r['block_edge'],
			 r['raw_reads_10kb'], round(r['mcov'], 2), int(r['linked']), r['neigh'], r['neigh10kb'],
			 r['ref_bases'], r["plus"], r["minus"], round(r["strand_binom_t"], 4), r['n_gaps'], round(r["n_sa"], 2),
			 round(r["n_xa"], 2),
			 round(r["n_unmapped_mates"], 2), r["double_clips"], r["remap_score"], r["remap_ed"], r["bad_clip_count"],
			 round(r["fcc"], 3), r["n_small_tlen"], r["ras"], r['fas'],
			 round(r["inner_cn"], 3), round(r["outer_cn"], 3), round(r["compress"], 2), round(r["ref_rep"], 3), round(r["jitter"], 3), r['prob']
			 ]
		return v

	else:
		v = [r["GT"], r["GQ"], r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
			 r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
			 r['sc'], r['bnd'], round(r['sqc'], 2), round(r['scw'], 1), round(r['clip_qual_ratio'], 3), r['block_edge'], r['raw_reads_10kb'], round(r['mcov'], 2), int(r['linked']), r['neigh'], r['neigh10kb'],
			 r['ref_bases'], r["plus"], r["minus"], round(r["strand_binom_t"], 4), r['n_gaps'], round(r["n_sa"], 2),
			 round(r["n_xa"], 2), round(r["n_unmapped_mates"], 2), r["double_clips"], r["remap_score"], r["remap_ed"], r["bad_clip_count"], round(r["fcc"], 3), r["n_small_tlen"], r["ras"], r['fas'],
			 round(r["inner_cn"], 3), round(r["outer_cn"], 3), round(r["compress"], 2), round(r["ref_rep"], 3), round(r["jitter"], 3), r['prob']
			 ]
		return v


def gen_format_fields(r, df, names, extended, n_fields, small_output):

	if len(names) == 1:
		return {0: get_fmt(r, extended, small_output)}, {}

	cols = {}
	if "partners" in r:
		if not isinstance(r["partners"], set):
			if len(r["partners"]) == 0 or pd.isna(r["partners"]):
				r["partners"] = []
			else:
				r["partners"] = [int(i.split(",")[1]) for i in r["partners"].split("|")]
		for idx in r["partners"]:
			if idx in df.index:  # might be already dropped
				r2 = df.loc[idx]
				cols[r2["table_name"]] = r2

	if "table_name" in r:
		cols[r["table_name"]] = r
	
	format_fields = sortedcontainers.SortedDict()
	for name in names:
		if name in cols:
			row = cols[name]
			format_fields[name] = get_fmt(row, extended, small_output)
		else:
			format_fields[name] = [0] * n_fields

	return format_fields, cols



def to_vcf(df, args, names, outfile, show_names=True,  contig_names="", extended_tags=False, header=None,
		   small_output_f=True, sort_output=True, n_fields=None):

	if header is None:
		HEADER = get_headers(extended_tags)
	else:
		HEADER = header

	outfile.write(HEADER+"\t"+names[0]+"\n")#.format(contig_names) + "\t" + "\t".join(names) + "\n")

	if show_names:
		logging.info("Input samples: {}".format(str(list(names))))

	version = pkg_resources.require("dysgu")[0].version

	seen_idx = set([])
	cnames = ['raw_reads_10kb', 'NMpri', 'NMsupp', 'MAPQpri', 'MAPQsupp', "NMbase", "n_gaps"]
	if extended_tags:
		cnames += ['DP', 'DN', 'DApri', 'DAsupp']

	for col in cnames:
		if col in df.columns:
			df[col] = df[col].astype(float).round(2)

	for col in ['maxASsupp', 'neigh', 'neigh10kb']:
		if col in df.columns:
			df[col] = [int(i) if (i == i and i is not None) else 0 for i in df[col]]

	count = 0
	recs = []

	if args is not None:
		add_kind = args["add_kind"] == "True"
		if args["metrics"]:
			small_output_f = False
		if args["verbosity"] == '0':
			df["contigA"] = [''] * len(df)
			df["contigB"] = [''] * len(df)
		elif args["verbosity"] == '1':
			has_alt = [True if isinstance(i, str) and i[0] != '<' else False for i in df['variant_seq']]
			df["contigA"] = ['' if a else c for c, a in zip(df['contigA'], has_alt)]
			df["contigB"] = ['' if a else c for c, a in zip(df['contigB'], has_alt)]

		n_fields = len(col_names(extended_tags, small_output_f)[-1])
	else:
		small_output_f = False

	for idx, r in df.iterrows():
		if idx in seen_idx:
			continue
		#print(r, df, names, extended_tags, n_fields, small_output_f)
		format_f, df_rows = gen_format_fields(r, df, names, extended_tags, n_fields, small_output_f)
		#print(format_f)
		if "partners" in r and r["partners"] is not None and r["partners"] != ".":
			seen_idx |= set(r["partners"])

		r_main = make_main_record(r, version, count, format_f, df_rows, add_kind, extended_tags, small_output_f)
		#print(r_main)
		recs.append(r_main)
		count += 1

	if sort_output:
		for rec in sorted(recs, key=lambda x: (x[0], x[1])):
			outfile.write("\t".join(list(map(str, rec))) + "\n")
	else:
		for rec in recs:
			outfile.write("\t".join(list(map(str, rec))) + "\n")

	return count



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


		if os.path.isfile(filename):
			if filename.lower().endswith((".vcf", ".csv")):
				files = [filename]
			else:
				raise ValueError(f"{filename[-3:]} invalid file type")
		if os.path.isdir(filename):
			files = os.listdir(filename)
			print("all", files)
			if os.path.isdir(outdir):
				already_done = os.listdir(outdir)
				print(already_done)
				if len(already_done) > 0:
					already_done = [i.replace("_filtered", "") for i in already_done]
					print("already", already_done)
					files = list(set(files) - set(already_done))
					print("now doing", files)
			files = [os.path.join(filename, i) for i in files if i.lower().endswith((".vcf", ".csv"))]
			if len(files) == 0:
				print("No vcf or csv files")
				exit

		merged_fmt = merged[-3:]
		if merged_fmt == "csv":
			print("Reading csv and converting dtypes")
			merged_df = pd.read_csv(merged)
			col_map = {"chrA": ("chrA", str),
			   "posA": ("posA", np.int64),
			   "chrB": ("chrB", str),
			   "grp_id": ("grp_id", np.int64),
			   "n_in_grp": ("n_in_grp", np.int64),
			   "posB": ("posB", np.int64),
			   "posB_tra": ("posB_tra", np.int64),
			   "sample": ("sample", str),
			   "event_id": ("event_id", np.int64),
			   "kind": ("kind", str),
			   "svtype": ("svtype", str),
			   "join_type": ("join_type", str),
			   "cipos95A": ("cipos95A", np.int64),
			   "cipos95B": ("cipos95B", np.int64),
			   "NMpri": ("NMpri", float),
			   "NMbase": ("NMbase", float),
			   "NMsupp": ("NMsupp", float),
			   "MAPQpri": ("MAPQpri", float),
			   "MAPQsupp": ("MAPQsupp", float),
			   "NP": ("NP", np.int64),
			   "query_overlap": ("query_overlap", np.int64),
			   "maxASsupp": ("maxASsupp", np.int64),
			   "su": ("su", np.int64),
			   "spanning": ("spanning", np.int64),
			   "pe": ("pe", np.int64),
			   "supp": ("supp", np.int64),
			   "sc": ("sc", np.int64),
			   "bnd": ("bnd", np.int64),
			   "sqc": ("sqc", float),
			   "scw": ("scw", float),
			   "clip_qual_ratio": ("clip_qual_ratio", float),
			   "type": ("type", str),
			   "block_edge": ("block_edge", np.int64),
			   "raw_reads_10kb": ("raw_reads_10kb", float),
			   "mcov": ("mcov", float),
			   "linked": ("linked", np.int64),
			   "contigA": ("contigA", str),
			   "contigB": ("contigB", str),
			   "ref_seq": ("ref_seq", str),
			   "variant_seq": ("variant_seq", str),
			   "gc": ("gc", float),
			   "neigh": ("neigh", np.int64),
			   "neigh10kb": ("neigh10kb", np.int64),
			   "rep": ("rep", float),
			   "rep_sc": ("rep_sc", float),
			   "svlen_precise": ("svlen_precise", np.int64),
			   "n_expansion": ("n_expansion", np.int64),
			   "stride": ("stride", np.int64),
			   "exp_seq": ("exp_seq", str),
			   "ref_poly_bases": ("ref_poly_bases", np.int64),
			   "GT": ("GT", str),
			   "GQ": ("GQ", int),
			   "ref_bases": ("ref_bases", np.int64),
			   "svlen": ("svlen", np.int64),
			   "plus": ("plus", np.int64),
			   "minus": ("minus", np.int64),
			   "strand_binom_t": ("strand_binom_t", float),
			   "prob": ("prob", float),
			   "n_gaps": ("n_gaps", float),
			   "n_sa": ("n_sa", float),
			   "n_xa": ("n_xa", float),
			   "n_unmapped_mates": ("n_unmapped_mates", np.int64),
			   "double_clips": ("double_clips", np.int64),
			   "remap_score": ("remap_score", np.int64),
			   "remap_ed": ("remap_ed", np.int64),
			   "bad_clip_count": ("bad_clip_count", np.int64),
			   "n_small_tlen": ("n_small_tlen", np.int64),
			   "ras": ("ras", np.int64),
			   "fas": ("fas", np.int64),
			   "inner_cn": ("inner_cn", float),
			   "outer_cn": ("outer_cn", float),
			   "compress": ("compress", float),
			   "fcc": ("fcc", float),
			   "ref_rep": ("ref_rep", float),
			   "jitter": ("jitter", float),
			   "left_ins_seq": ("left_ins_seq", str),
			   "right_ins_seq": ("right_ins_seq", str),
			   }
			merged_df.rename(columns={k: v[0] for k, v in col_map.items()}, inplace=True)
			merged_df["GQ"] = pd.to_numeric(merged_df["GQ"], errors='coerce').fillna(".")
			merged_df["GQ"] = merged_df["GQ"].replace(".", 0, regex=True)
			merged_df["GQ"] = merged_df["GQ"].replace(np.nan, 0, regex=True)
			for k, dtype in col_map.values():
				if k in merged_df:
					if merged_df[k].dtype != dtype:
						try:
							merged_df[k] = merged_df[k].astype(dtype)
						except ValueError or OverflowError:
							raise ValueError("Problem for feature {}, could not intepret as {}".format(k, dtype))

						if dtype == str:
							merged_df[k] = merged_df[k].fillna("")
						else:
							merged_df[k] = merged_df[k].fillna(0)
						#try:
						#	merged_df[k] = merged_df[k].astype(dtype)
						#except ValueError or OverflowError:
						#	#print(list(df[k][0]))
						#	raise ValueError("Problem for feature {}, could not intepret as {}".format(k, dtype))
			#if "contigA" not in merged_df:
			#	merged_df["contigA"] = [""] * len(merged_df)
			#if "contigB" not in merged_df:
			#	merged_df["contigB"] = [""] * len(merged_df)

			#if 'posB_tra' in merged_df:
			#	merged_df["posB"] = [i if svt != 'TRA' else j for i, j, svt in zip(merged_df['posB'], merged_df['posB_tra'], merged_df['svtype'])]
			#	del merged_df['posB_tra']
			merged_df["posA"] = merged_df["posA"].astype(int) # convert to 0-indexing
			merged_df["posB"] = merged_df["posB"].astype(int)
			merged_df = merged_df.drop(columns=["partners"])
			print(f"Merged file has {len(merged_df)} variants occuring accross all samples")
			
		elif merged_fmt == "vcf":
			approx_vars = file_len(merged) - check_header_len(merged) - 1
			print(f"Converting vcf to csv (expect {approx_vars} variants)")
			merged_df, merged_header, merged_fields, merged_num_vars = combined_vcf_to_df(merged)
			print(f"Merged file has {merged_num_vars} variants (occuring {len(merged_df)} times)")
		else:
			raise ValueError(f"{merged_fmt} invalid file type")
		#merged_df["table_name"] = "common"

				
		for f in files:
			sample = str(os.path.basename(f).rsplit(".", 1)[0]) 
			out_file_path = os.path.join(outdir, sample+suffix+"."+fmt)

			sample_fmt = f[-3:]
			if sample_fmt == "csv":
				in_df = pd.DataFrame.from_csv(f)
			elif sample_fmt == "vcf":
				in_df, in_header, in_fields = vcf_to_df(f)
			else:
				print(f"{sample_fmt} invalid file type")	
			before = len(in_df)
			print(f"Input file {sample} has {len(in_df)} variants")
			#in_df["table_name"] = sample

			print(f"Attempting merge {sample}...")
			t0 = time.time()
			all_df = dysgu.merge_dysgu_df(in_df, merged_df)
			all_df = all_df.loc[:,~all_df.columns.duplicated()].copy()
			print(f"Merge took {' '.join([i+z for i, z in zip(str(timedelta(seconds=time.time()-t0)).split(':'), ['h', 'm', 's'])])}")
			#print(all_df)
			all_df.to_csv("all.csv")
			filtered_df = all_df.loc[(all_df["partners"].map(len) == 0) & (all_df["sample"] != "common")]
			filtered_df["table_name"] = sample
			if fmt == "vcf":
				with open(out_file_path, "w") as of:
					to_vcf(filtered_df, {"add_kind": True, "metrics": True, "verbosity": "1"}, [sample], of, header=in_header)
			if fmt == "csv":
				filtered_df.to_csv(out_file_path)
			print(f"Filtered file has {len(filtered_df)} variants (removed {before - len(filtered_df)})")


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
