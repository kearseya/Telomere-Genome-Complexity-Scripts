import os
import pandas as pd

samples = os.listdir("icgc_download_data")
for d in samples:
	fn = os.path.join("icgc_download_data", d, "copy_number_somatic_mutation.tsv")
	print(fn)
	if os.path.exists(fn):
		data = pd.read_csv(fn, sep="\t")
		print(set(list(data["icgc_specimen_id"])))
