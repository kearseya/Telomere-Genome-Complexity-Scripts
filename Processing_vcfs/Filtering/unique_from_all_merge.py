#!/usr/bin/python3
import pandas as pd
import os
import gzip
import sys
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

@click.command()
@click.argument("path")
@click.argument("sample")
@click.option("-o", default=".")
def filter(path, sample, o):
	fin = read_from_inputfile(path)
	# Parse sample
	header = ""
	last = ""
	keep_individual_meta = False
	header_length = 0
	for line in fin:
		if line[:2] == "##":
			header += line
			header_length += 1
		else:
			break
			#header += "\t".join(line.split("\t")[:9])
			#last = line
			#break

	samples = last.strip().split("\t")[9:]
	n_samples = len(samples)

	if path == '-':
		path_h = sys.stdin
	else:
		path_h = path

	df = pd.read_csv(path_h, index_col=None, sep="\t", header=header_length)
	combined_count = len(df)
	df = df[df[sample] != "0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0"]
	sample_count = len(df)
	c = list(df.columns)[9:]
	c.remove(sample)
	ids = set()
	for s in c:
		indexes = df.index[df[s] != "0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0"].tolist()
		for i in indexes:
			ids.add(i)
	df = df.drop(ids, axis=0)
	df = df.drop(c, axis=1)
	filtered_count = len(df)
	print(f"{sample} contains {sample_count}/{combined_count} variants of which {filtered_count} ({len(df[df['FILTER'] == 'PASS'])}/{len(df[df['FILTER'] == 'lowProb'])}) are unique (removed {sample_count-filtered_count})")
	outname = os.path.join(o, f"{sample}_unique.vcf")
	df.to_csv(outname, sep="\t", index=False)

if __name__ == "__main__":
	filter()
