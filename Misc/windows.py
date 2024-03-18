import pandas as pd
import sys
from rich.progress import Progress
import fileinput

df = pd.read_csv("hg38_prep_10k.bed", sep="\t", header=None)
df.columns = ["chrom", "start", "end"]
df["depth"] = 0

sample = sys.argv[1]

d = {}
chroms = set()

for i in df.groupby("chrom"):
    d[i[0]] = []
    chroms.add(i[0])
    for r in i[1].itertuples(index=False, name=None):
        d[i[0]].append(list(r)[1:])

total_spaces = 0
for c in d:
    total_spaces += int(d[c][-1][1])
print(total_spaces)

total_spaces = total_spaces/100

with Progress() as progress:
    task = progress.add_task("Getting coverage: ", total=total_spaces)
    for line in fileinput.input("-", encoding="utf-8"):
        l = line.split()
        if len(l) > 0:
            if l[0] in chroms:
                d[l[0]][int(l[1])//10000][-1] += int(l[2])
                progress.update(task, advance=0.01)


for c in d:
    for i, _ in enumerate(d[c]):
        d[c][i][-1] = d[c][i][-1]/(int(d[c][i][1])-int(d[c][i][0]))

rec = []
for c in d:
    for i in d[c]:
        rec.append([c]+i)

df = pd.DataFrame.from_records(rec)
df.to_csv(f"coverage/{sample}_cov.bed", sep="\t", index=False, quoting=None)
print(df)

