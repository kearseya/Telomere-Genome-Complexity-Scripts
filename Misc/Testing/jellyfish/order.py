import array
import re

#for line in open("test_mapping"):
#   print(line)

#order = array.array("B", [])
a = array.array("I", [])
for i in range(1, 23):
    pos = 0
    order = []
    count = 0
    #print("")
    #print(i)
    for line in open("hg19_mappings"):
        #print(pos, end="\r")
        l = re.findall('\d+', line)
        if len(l) == 0:
            if line == "*":
                pos += 1
            continue
        if int(l[0]) == i:
            if "alt" not in line:
                if len(l) > 2:
                    pos += 1
                    continue
                if len(l) == 2:
                    l = (int(l[1]), pos)
                    order.append(l)
                    count += 1
        pos += 1
    
    inorder = sorted(order, key=lambda x: x[0])
    print(inorder[:25])
    for z in inorder:
        a.append(z[1])
    print(i, count)


for i in range(23, 25):
    pos = 0
    order = []
    count = 0
    #print(i)
    for line in open("hg19_mappings"):
        #print(pos, end="\r")
        l = re.findall('\d+', line)
        if i in {23, 24}:
            if len(l) == 1:
                if line[3] == "X":
                    l.insert(0, 23)
                if line[3] == "Y":
                    l.insert(0, 24)
        if len(l) == 0:
            if line == "*":
                pos += 1
            continue
        if int(l[0]) == i:
            if "alt" not in line:
                if len(l) > 2:
                    pos += 1
                    continue
                if len(l) == 2:
                    #l[1] = int(l[1])
                #l.append(pos)
                    l = (int(l[1]), pos)
                    order.append(l)
                    count += 1
            #else:
                #order.append((1, pos))
        pos += 1
    
    inorder = sorted(order, key=lambda x: x[0])
    #print(inorder[:25])
    for z in inorder:
        a.append(z[1])
    print(i, count)


a.tofile(open("hg19_index", "wb"))

