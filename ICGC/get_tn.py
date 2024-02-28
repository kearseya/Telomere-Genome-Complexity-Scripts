import pandas as pd

d = pd.read_csv("breast_WGS.tsv", sep="\t")
x = d[["File Name", "Specimen Type"]]

## rename columns
x = x.rename(columns = {"File Name": "file", "Specimen Type": "type"})

## filter mini and pcawg
x = x[x["file"].str.contains("mini|PCAWG")==False]
x = x.reset_index()

## remove extension
x["file"] = x.file.apply(lambda i: i.split('.')[0])


## change type to one word
x.loc[x["type"].str.contains("Normal"), "type"] = "normal"
x.loc[x["type"].str.contains("tumour"), "type"] = "tumour"

## write file
x.to_csv("icgc_tn.csv", index=False)
