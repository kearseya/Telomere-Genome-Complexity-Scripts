import os
import matplotlib.pyplot as plt
import pandas as pd

for i in os.listdir("beds"):
    if i.endswith(".bed") == False:
        continue
    df = pd.read_csv(os.path.join("beds", i), header=None, sep="\t")
    df = df.rename(columns={1: "start", 2: "end"})
    df["len"] = df["end"] - df["start"]
    print(i)
    print(df[df["len"] < 0])
    df = df.drop(df.index[df["len"] < 0])
    plt.boxplot(df["len"])
    #plt.show()
plt.show()
