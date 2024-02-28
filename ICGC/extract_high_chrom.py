import re

labels = {"Cancer type", "Position", "Type", "Interleaved intrachr. SVs", "Total SVs (intrachr. + transl.)", "SV types", "SVs in sample", 
        "Oscillating CN (2 and 3 states)", "CN segments", "FDR fragment joints", "FDR chr. breakp. enrich.", "Linked to chrs", "Purity, ploidy"}

l = []

# add raw lines
with open("high_chrom.txt", "r") as f:
    for line in f:
        line = line.strip()
        if line not in labels:
            l.append(line)

#print(l)

# remove left columns
for p, i in enumerate(l):
    for n in labels:
        if n in str(i):
            l[p] = i.split(n)[1].strip()

# seperate name from purity and plurity
name_purity_split = re.compile(r"([0-9]\.[0-9]*\, [0-9]\.[0-9][0-9])+").split
c = l.copy()
l = [s for i in c for s in name_purity_split(i) if s]
name_na_split = re.compile(r"(NA, NA)+").split
c = l.copy()
l = [s for i in c for s in name_purity_split(i) if s]
print(l)

d = {}
pos = 1

#sv_type = re.findall(r"([0-9];)")

# table structure: total 19
# name, cancer type, position, type, interleaved intrachr sv, total sv, #
# sv del, sv dup, sv h2hlNV. sv t2tlNV, sv tra
# svs, osc cn 2, osc cn 3, cn segments, fdr frag join, fdr chr breakp, purity, ploidy
for i in l:
    if pos == 1:
        d[i] = {}
        s = i
    if pos == 2:
        d[s]["cancer type"] = i
    if pos == 3:
        d[s]["position"] = i
    if pos == 4:
        d[s]["type"] = i
    if pos == 5:
        d[s]["interleaved"] = i
    if pos == 6:
        d[s]["total sv"] = i
    if pos == 7:
        if "DEL" in i:
            for x, z in enumerate(re.findall(r"([0-9]*;)", i)):
                if x == 0:
                    d[s]["sv del"] = z.split(";")[0]
                if x == 1:
                    d[s]["sv dup"] = z.split(";")[0]
                if x == 2:
                    d[s]["sv h2hlNV"] = z.split(";")[0]
            pos = pos-1
        if "t2tINV" in i:
            for x, z in enumerate(re.findall(r"(\s[0-9]*)", i)):
                if x == 0:
                    d[s]["sv t2tINV"] = z
                if x == 2:
                    d[s]["sv tra"] = z
    if pos == 8:
        d[s]["SVs"] = i
    if pos == 9:
        if "," in i:
            d[s]["osc cn 2"] = i.split(",")[0]
            d[s]["osc cn 3"] = i.split(",")[1].strip()
    if pos == 10:
        d[s]["cn seg"] = i
    if pos == 11:
        d[s]["fdr frag"] = i
    if pos == 12:
        d[s]["fdr break"] = i
        d[s]["links"] = []
    if pos >= 13:
        if ":" in i:
            d[s]["links"].append(i)
        if "," in i:
            d[s]["purity"] = i.split(",")[0]
            d[s]["ploidy"] = i.split(",")[1].strip()
            pos = 0

    pos += 1


print(d)





