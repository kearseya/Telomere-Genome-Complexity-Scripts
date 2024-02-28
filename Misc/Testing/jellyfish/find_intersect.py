hg38 = open("38", "r")
out = open("kmers", "w")

m = 0
c = 0
for p, i in enumerate(hg38):
    hg19 = open("19", "r")
    for _ in range(p+c):
        next(hg19)
    for z, x in enumerate(hg19):
        if i == x:
            m += 1
            c += z
            out.write(i)
            break
        if z > p+2500:
            break

# max = 3286714

hg19.close()
hg38.close()
out.close()
