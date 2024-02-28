import pandas as pd

df = pd.read_csv("all_files.tsv", sep="\t")
df = df[df["Format"] == "BAM"]
df = df.reset_index()

l = []
with open("nomini_object_ids", "r") as f:
    for line in f:
        l.append(line.strip())

df = df[df["Object ID"].isin(l)]
df = df.reset_index()

df = df[["File ID", "Object ID", "File Name", "ICGC Donor", "Specimen ID", "Specimen Type", "Sample ID"]]

norm = df[df["Specimen Type"].str.startswith("Normal")]
norm = norm[["File Name", "Object ID", "ICGC Donor"]]
norm = norm.rename(columns={"File Name": "normal", "Object ID": "nID"})
tumr = df[df["Specimen Type"].str.startswith("Primary")]
tumr = tumr[["File Name", "Object ID", "ICGC Donor"]]
tumr = tumr.rename(columns={"File Name": "tumour", "Object ID": "tID"})

merge = pd.merge(norm, tumr, on = "ICGC Donor")
merge.to_csv("icgc_tn_on_donor.csv")



