import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

names = pd.read_csv("breast_cancer_ICGC/breast_WGS.tsv", sep="\t")
names.columns = names.columns.str.replace(' ', '_', regex=True)
names = names[~names.File_Name.str.contains("mini")]
names = names[~names.File_Name.str.contains("PCAWG")]
samples = os.listdir("icgc_download_data")

fntsn = pd.Series(names["ICGC_Donor"].values, index=names["File_Name"])#.to_dict()
#fntsn = {k:v for k,v in fntsn.items() if not 'mini' in k}
#fntsn = {k:v for k,v in fntsn.items() if not 'PCAWG' in k}
fntsn = {k.split(".bam")[0]:v for k,v in fntsn.items()}

dids = list(set(list(names["ICGC_Donor"])))
#with pd.option_context('display.max_rows', None):#, 'display.max_columns', None):
#	print(names[["File_Name", "Specimen_ID", "ICGC_Donor"]])
sntsi = {}
sntfn = {}
for i in dids:
	sntsi[i] = {}
	sids = names[names["ICGC_Donor"] == i]
	sntsi[i]["normal"] = sids[sids["Specimen_Type"].apply(lambda x: "Normal" in x)]["Specimen_ID"].iloc[0]
	sntsi[i]["tumour"] = sids[sids["Specimen_Type"].apply(lambda x: "tumour" in x)]["Specimen_ID"].iloc[0]
	sntfn[i] = {}
	sntfn[i]["normal"] = sids[sids["Specimen_Type"].apply(lambda x: "Normal" in x)]["File_Name"].iloc[0].split(".bam")[0]
	sntfn[i]["tumour"] = sids[sids["Specimen_Type"].apply(lambda x: "tumour" in x)]["File_Name"].iloc[0].split(".bam")[0]

for d in [sntsi, sntfn]:
	for v in d:
		if len(d[v]) != 2:
			print("error")

print(len(sntsi))
print(len(sntfn))

pred = pd.read_csv("/path/to/project/TL_prediction/allpredot.csv")
pred = pred[["sample", "g", "short", "long_conf", "short_conf"]]
pred = pred[pred["g"] == "i"]

s = list(sntsi.keys())
print(sorted(s))
d = {}
df = {}
#print(fntsn["6a14cdd22dacb3527dd3664837577ff3"])
for n, i in enumerate(s):
	print(f"{n}/{len(s)}", end="\r")
	if i in samples:
		d[i] = {}
		#tmptype = names[names["ICGC_Donor"] == i]["Specimen_Type"].iloc[0]
		#d[i]["type"] = tmptype.split(" ")[0]
		d[i]["normal"] = {}
		d[i]["tumour"] = {}
		d[i]["normal"]["spid"] = sntsi[i]["normal"]
		d[i]["tumour"]["spid"] = sntsi[i]["tumour"]
		d[i]["normal"]["sample"] = sntfn[i]["normal"]
		d[i]["tumour"]["sample"] = sntfn[i]["tumour"]
		files = os.listdir(os.path.join("icgc_download_data", i))
		for fn in files:
			if os.path.basename(fn) == "copy_number_somatic_mutation.tsv":
				f = pd.read_csv(os.path.join("icgc_download_data", i, "copy_number_somatic_mutation.tsv"), sep="\t")
				d[i]["normal"]["cnsm"] = len(f[f["icgc_specimen_id"] == sntsi[i]["normal"]].index)
				d[i]["tumour"]["cnsm"] = len(f[f["icgc_specimen_id"] == sntsi[i]["tumour"]].index)
			if os.path.exists(os.path.join("icgc_download_data", i, "copy_number_somatic_mutation.tsv")) == False:
				d[i]["normal"]["cnsm"] = np.NaN
				d[i]["normal"]["cnsm"] = np.NaN
			if os.path.basename(fn) == "simple_somatic_mutation.open.tsv":
				f = pd.read_csv(os.path.join("icgc_download_data", i, "simple_somatic_mutation.open.tsv"), sep="\t")
				d[i]["normal"]["ssm"] = len(f[f["icgc_specimen_id"] == sntsi[i]["normal"]].index)
				d[i]["tumour"]["ssm"] = len(f[f["icgc_specimen_id"] == sntsi[i]["tumour"]].index)
			if os.path.exists(os.path.join("icgc_download_data", i, "simple_somatic_mutation.open.tsv")) == False:
				d[i]["normal"]["ssm"] = np.NaN
				d[i]["tumour"]["ssm"] = np.NaN

for i in d:
	print(i)
	for t in ["normal", "tumour"]:
		fn = d[i][t]["sample"]
		df[fn] = {}
		df[fn]["did"] = i
		df[fn]["type"] = t
		for v in ["spid", "cnsm", "ssm"]:
			if v in d[i][t]:
				df[fn][v] = d[i][t][v]
			else:
				df[fn][v] = np.NaN

print(df)

d = pd.DataFrame(df).transpose()
d["sample"] = d.index
d = d.reset_index()
#d = d[["sample", "type", "cnsm", "ssm"]]
with pd.option_context('display.max_rows', None):#, 'display.max_columns', None):
	print(d)
#print(fntsn)
pred = pd.read_csv("/path/to/project/TL_prediction/allpredot.csv")
pred = pred[["sample", "g", "short", "long_conf", "short_conf"]]

d = pd.merge(pred, d, how="left", on="sample")

with pd.option_context('display.max_rows', None):#, 'display.max_columns', None):
	print(d)

d.to_csv("mutation_counts.csv")

#print(d)
short_cnsm = d[d["short"] == True]["cnsm"].dropna()
long_cnsm = d[d["short"] == False]["cnsm"].dropna()
short_ssm = d[d["short"] == True]["ssm"].dropna()
long_ssm = d[d["short"] == False]["ssm"].dropna()

norm_cnsm = d[d["type"] == "normal"]["cnsm"].dropna()
tum_cnsm = d[d["type"] == "tumour"]["cnsm"].dropna()
tum_cnsm = tum_cnsm.drop(pd.to_numeric(tum_cnsm).idxmin())
norm_ssm = d[d["type"] == "normal"]["ssm"].dropna()
tum_ssm = d[d["type"] == "tumour"]["ssm"].dropna()

#bins=np.histogram(np.hstack((short_cnsm, long_cnsm)), bins=60)[1]
#plt.hist(short_cnsm, bins, alpha=0.5)
#plt.hist(long_cnsm, bins, alpha=0.5)
#plt.legend(["short", "long"])
#plt.xlabel("# copy number somatic mutations")
#plt.show()
#
#bins=np.histogram(np.hstack((short_ssm, long_ssm)), bins=100)[1]
#plt.hist(short_ssm, bins, alpha=0.5)
#plt.hist(long_ssm, bins, alpha=0.5)
#plt.legend(["short", "long"])
#plt.xlabel("# simple somatic mutations")
#plt.show()
#
#bins=np.histogram(np.hstack((norm_cnsm, tum_cnsm)), bins=60)[1]
#plt.hist(norm_cnsm, bins, alpha=0.5)
#plt.hist(tum_cnsm, bins, alpha=0.5)
#plt.legend(["norm", "tum"])
#plt.xlabel("# copy number somatic mutations")
#plt.show()
#
#bins=np.histogram(np.hstack((norm_ssm, tum_ssm)), bins=100)[1]
#plt.hist(norm_ssm, bins, alpha=0.5)
#plt.hist(tum_ssm, bins, alpha=0.5)
#plt.legend(["norm", "tum"])
#plt.xlabel("# simple somatic mutations")
#plt.show()


shorts = d[d["short"] == True]
longs = d[d["short"] == False]
normals = d[d["type"] == "normal"]
tumours = d[d["type"] == "tumour"]

#print(shorts)
#print(longs)
#print(normals)
#print(tumours)


cnsm_data = {}
for sl, i in zip(["short", "long"], [shorts, longs]):
	for nt, x in zip(["normal", "tumour"], [normals, tumours]):
		tmp = pd.merge(i, x[["sample"]], how="inner", on="sample")
		tmp = tmp.loc[:,~tmp.columns.duplicated()]
		#print(tmp)
		cnsm_data[f"{sl}_{nt}"] = tmp["cnsm"].dropna()

ssm_data = {}
for sl, i in zip(["short", "long"], [shorts, longs]):
	for nt, x in zip(["normal", "tumour"], [normals, tumours]):
		tmp = pd.merge(i, x[["sample"]], how="inner", on="sample")
		#tmp = tmp.loc[:,~tmp.columns.duplicated()]
		ssm_data[f"{sl}_{nt}"] = tmp["ssm"].dropna()


cnsm_data["short_tumour"] = cnsm_data["short_tumour"].drop(pd.to_numeric(cnsm_data["short_tumour"]).idxmin())

#print(data)
bins=np.histogram(np.hstack((cnsm_data["long_tumour"], cnsm_data["short_tumour"])), bins=100)[1]		#(cnsm_data["long_normal"], cnsm_data["short_normal"], 
							 
#plt.hist(cnsm_data["long_normal"], bins, alpha=0.7)
#plt.hist(cnsm_data["short_normal"], bins, alpha=0.7)
plt.hist(cnsm_data["long_tumour"], bins, alpha=0.7)
plt.hist(cnsm_data["short_tumour"], bins, alpha=0.7)
plt.legend(["long_tumour", "short_tumour"])#["long_normal", "short_normal", "long_tumour", "short_tumour"])
plt.xlabel("#cnsm")
plt.show()


bins=np.histogram(np.hstack((ssm_data["long_normal"], ssm_data["short_normal"], 
							 ssm_data["long_tumour"], ssm_data["short_tumour"])), bins=100)[1]
plt.hist(ssm_data["long_normal"], bins, alpha=0.7)
plt.hist(ssm_data["short_normal"], bins, alpha=0.7)
plt.hist(ssm_data["long_tumour"], bins, alpha=0.7)
plt.hist(ssm_data["short_tumour"], bins, alpha=0.7)
plt.legend(["long_normal", "short_normal", "long_tumour", "short_tumour"])
plt.xlabel("#ssm")
plt.show()


	

#for n, i in enumerate(s):
#	print(f"{n}/{len(s)}", end="\r")
#	if fntsn[i] in samples:
#		d[sntfn[i]] = {}
#		tmptype = names[names["ICGC_Donor"] == i]["Specimen_Type"].iloc[0]
#		d[sntfn[i]]["type"] = tmptype.split(" ")[0]
#		files = os.listdir(os.path.join("icgc_download_data", i))
#		for f in files:
#			if os.path.basename(f) == "copy_number_somatic_mutation.tsv":
#				d[sntfn[i]]["cnsm"] = len(pd.read_csv(os.path.join("icgc_download_data", i, "copy_number_somatic_mutation.tsv")).index)
#			if os.path.exists(os.path.join("icgc_download_data", i, "copy_number_somatic_mutation.tsv")) == False:
#				d[sntfn[i]]["cnsm"] = np.NaN
#			if os.path.basename(f) == "simple_somatic_mutation.open.tsv":
#				d[sntfn[i]]["ssm"] = len(pd.read_csv(os.path.join("icgc_download_data", i, "simple_somatic_mutation.open.tsv")).index)
#			if os.path.exists(os.path.join("icgc_download_data", i, "simple_somatic_mutation.open.tsv")) == False:
#				d[sntfn[i]]["ssm"] = np.NaN
#
#		#try:
#		#	d[sntfn[i]]["cnsm"] = len(pd.read_csv(os.path.join("icgc_download_data", i, "copy_number_somatic_mutation.tsv")).index)
#		#except:
#		#	continue
#		#try:
#		#	d[sntfn[i]]["ssm"] = len(pd.read_csv(os.path.join("icgc_download_data", i, "simple_somatic_mutation.open.tsv")).index)
#		#except:
#		#	continue
#

