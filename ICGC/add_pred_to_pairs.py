import pandas as pd

p = pd.read_csv("test.csv")
l = pd.read_csv("tmp_icgc_lin_pred.csv")
p["tumour"] = p["tumour"].str.replace(".bam", "")
print(p["tumour"])
print(l["sample"])

l = l[["sample", "pred_val"]]
l = l.rename({"sample": "tumour"}, axis=1)

df = p.merge(l, how="left", on="tumour")
df = df.dropna()
print(df)
df.to_csv("icgc_pairs.csv")
