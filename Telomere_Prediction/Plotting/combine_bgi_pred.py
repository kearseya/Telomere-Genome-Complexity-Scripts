import pandas as pd
import numpy as np

b = pd.read_csv("bgi_stela.csv")
bq = pd.read_csv("bgi_qmotif.csv")
bc = pd.read_csv("bgi_telomerecat_length.csv")
bt = pd.read_csv("bgi_tumour_predicted_tl.csv")
bs = pd.read_csv("bgi_telseq.csv")

b = b.merge(bq, on="sample", how="left")
bc = bc.rename({"Length": "telomerecat"}, axis=1)
b = b.merge(bc[["sample", "telomerecat"]], on="sample", how="left")
bt = bt.rename({"Tumour": "sample"}, axis=1)
b = b.merge(bt[["sample", "CompuTel"]], on="sample", how="left")
bs = bs.rename({"ReadGroup": "sample", "LENGTH_ESTIMATE": "telseq"}, axis=1)
b = b.merge(bs[["sample", "telseq"]], on="sample", how="left")

b["telomerecat"] = b["telomerecat"]/1000
b["qmotif"] = b["qmotif"]/1000 

print(b)
b.to_csv("bgi_all_tool_pred.csv", quoting=None, index=False)
