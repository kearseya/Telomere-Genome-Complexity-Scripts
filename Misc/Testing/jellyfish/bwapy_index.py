from bwapy import BwaAligner
import joblib
import array
import numpy as np
import re

a = BwaAligner("/scratch/ProjectDir/reference_genomes/hg38.fa")

count = 0
index_binary = array.array("B", [])
for i in open("/scratch/ProjectDir/jellyfish/filter.fa"):
	i = i.strip()
	if len(i) == 21:
			count += 1
			al = a.align_seq(i)
			if len(al) == 1:
				rname = al[0].rname
				l = re.findall('\d+', rname)
				if len(l) == 1:
					index_binary.append(1)
				if len(l) > 1:
					index_binary.append(0)
				if len(l) == 0:
					xy = rname[-1]
					if xy in {"X", "Y"}:
						index_binary.append(1)
					else:
						index_binary.append(0)
			else:
				index_binary.append(0)

	#print(count, end="\r")

index_binary.tofile(open("index_binary", "wb"))

