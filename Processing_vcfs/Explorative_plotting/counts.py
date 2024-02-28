#!/usr/bin/python3
import click
import matplotlib.pyplot as plt
#from collections import defaultdict
import collections
import operator
import re
import os
import pandas as pd

@click.command()
@click.argument("indir", nargs=1)
@click.argument("out", nargs=1)
def count(indir, out):
    if os.path.isfile(indir):
        if indir.lower().endswith((".vcf", ".csv")):
            files = [indir]
        else:
            print("Cannot process file format")
            exit
    elif os.path.isdir(indir):
        files = os.listdir(indir)
        files = [i for i in files if i.lower().endswith((".vcf", ".csv"))]
        files = [os.path.join(indir, i) for i in files]
    failed = []
    total = 0
    counts = {}
    types = ["DEL", "INS", "INV", "DUP", "TRA"]
    for i in files:
        total += 1
        fmt = i[-3:]
                
        if fmt == "vcf":
            counts[i] = {"lowProb": {}, "PASS": {}}
            for p in counts[i].keys():
                for t in types:
                    counts[i][p][t] = 0
            with open(i) as f:
                f.readline()
                for line in f:
                    if line.startswith("#"):
                        continue
                    line = line.rstrip().split()
                    #print(line[4], line[6])
                    try:
                        #print(i, line[6], str(re.search("SVTYPE=([A-Z]*)", line[7]).group(1)))
                        counts[i][line[6]][str(re.search("SVTYPE=([A-Z]*)", line[7]).group(1))] += 1
                    except:
                        pass
        
        elif fmt == "csv":
            counts[i] = {}
            df = pd.read_csv(i)
            # counts[i]["lowProb"] = dict(df[df["filter"] == "lowProb"]["svtype"].value_counts())
            # counts[i]["PASS"] = dict(df[df["filter"] == "PASS"]["svtype"].value_counts())
            counts[i]["PASS"] = dict(df["svtype"].value_counts())

        else:
            print(f"cannot process {i}, invalid format")

    #df = pd.DataFrame.from_dict(counts)
    df = pd.DataFrame.from_dict({(i,j): counts[i][j] 
                               for i in counts.keys() 
                               for j in counts[i].keys()},
                           orient='index')
    df.to_csv(out)
    df = pd.read_csv(out)
    df = df.rename(columns={"Unnamed: 0": "sample", "Unnamed: 1": "filter"})
    df["sample"] = [os.path.basename(i) for i in df["sample"]]
    df.to_csv(out, index=False)

if __name__ == "__main__":
    count()
#print(failed)
