import sys
import array
outf = sys.argv[1]
a = array.array("B", [])
for line in sys.stdin:
    v = int(line.strip().split(" ")[1])
    if v > 255:
        v = 255
    a.append(v)
a.tofile(open(outf, "wb"))

