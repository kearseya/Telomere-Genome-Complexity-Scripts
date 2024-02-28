#!/usr/bin/python3
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from sklearn.cluster import KMeans
import click


def get_len_array(df, vt, lpp):
    if lpp == "PASS":
        lpp = 1
    else:
        lpp = 0
    l = df[vt][0][lpp].split(",")
    print(l)
    l = list(filter(None, l))
    l[0] = l[0].replace("[", "")
    l[-1] = l[-1].replace("]", "")
    print(l)
    l = [int(i.strip(" '")) for i in l]
    return l


@click.command()
@click.argument("file")
@click.option("-l", "--length", default=None)
@click.option("-n", "--col-name", default="short")
@click.option("-s", "--sample-col", default="sample")
@click.option("-t", "--threshold", default=3.81)
@click.option("--scale", default="linear",
    type=click.Choice(["linear", "log", "symlog", "logit"]))
@click.option("-k", "--kgroups", default=5)
@click.option("-r", "--real-groups", default=None)
@click.option("--ksort", default=False, is_flag=True)
@click.option("--lsort", default=False, is_flag=True)
@click.option("--save", default=False, is_flag=True)
@click.option("--low-prob", default=False, is_flag=True)
@click.option("--bar", default=False, is_flag=True)
@click.option("--bar-lines", default=[None], multiple=True)
@click.option("--bar-numbers", default=[None], multiple=True)
@click.option("--bar-scale", default=True, is_flag=True)
@click.option("--box", default=True, is_flag=True)
@click.option("--hist", default=False, is_flag=True)
def plot_vars(file, length, col_name, sample_col, threshold, scale, 
                    kgroups, real_groups, ksort, lsort,
                    save, low_prob, 
                    bar, bar_lines, bar_numbers, bar_scale,
                    box, hist):
    df = pd.read_csv(file)
    # >>> import pandas as pd
    # >>> df = pd.read_csv("test.csv")
    # >>> l = df["DEL"][0].split(",")
    # >>> l[0] = l[0].replace("[", "")
    # >>> l[-1] = l[-1].replace("]", "")
    # >>> [int(i.strip(" '")) for i in l]

    var = df.columns[2:].tolist()
    print(var)
    df["total"] = df[var].sum(axis=1)
    df["sample"] = [i.replace("_filtered.vcf", "") for i in df["sample"].tolist()]
    df["sample"] = [i.replace("_unique.vcf", "") for i in df["sample"].tolist()]    
    df["sample"] = [i.replace("_markedSpec.vcf", "") for i in df["sample"].tolist()]    
    df["sample"] = [i.replace(".vcf", "") for i in df["sample"].tolist()]
    df["sample"] = [i.replace("_filtered.csv", "") for i in df["sample"].tolist()]
    df["sample"] = [i.replace("_unique.csv", "") for i in df["sample"].tolist()]    
    df["sample"] = [i.replace("_markedSpec.csv", "") for i in df["sample"].tolist()]    
    df["sample"] = [i.replace(".csv", "") for i in df["sample"].tolist()]

    ## filter by filter column ##
    if low_prob == False:
        df = df[df["filter"] == "PASS"]
    else:
        odf = df.copy()
        df = df.groupby(df["sample"], as_index=False).agg({"DEL": sum, "DEL": sum, "INS": sum, "INV": sum, "DUP": sum, "TRA": sum, "total": sum})#, col_name: "first"})
    ## reorder table
    if length == None:
        df = df.sort_values(by="total")
    else:
        ldf = pd.read_csv(length)
        if sample_col != "sample":
            ldf = ldf.rename({sample_col: "sample"})
        ldf = ldf[["sample", col_name]]
        df = pd.merge(df, ldf, on="sample")
        
        if col_name != "short" and df[col_name].dtype != bool:
            df["short"] = np.where(df[col_name] <= threshold, True, False)
        #print(df)
        if df[col_name].dtype == bool:
            df = df.sort_values(by=[col_name, "total"])
        else:
            if lsort == False:
                df = df.sort_values(by=["short", "total"])
            else:
                df = df.sort_values(by=col_name, ascending=False)
        if "length" in bar_lines:
            if df[col_name].dtype != bool:
                for j, z in enumerate(df[col_name]):
                    if z <= threshold:
                        lsplit = j
                        print(lsplit)
                        break
            else:
                for j, z in enumerate(df[col_name]):
                    if z == True:
                        lsplit = j
                        print("bool", lsplit)
                        break


    df = df.reset_index()
#   ## filter by filter column ##
#   if low_prob == False:
#       df = df[df["filter"] == "PASS"]
#   else:
#       odf = df.copy()
#       df = df.groupby(df["sample"], as_index=False).agg({"DEL": sum, "DEL": sum, "INS": sum, "INV": sum, "DUP": sum, "TRA": sum, "total": sum, col_name: "first"})


    ## add kmeans labels ##
    # kmeans = KMeans(n_clusters=kgroups, random_state=0).fit(df[var])
    # df["kmeans"] = kmeans.labels_
    #print(kmeans.labels_)
    #print(df["kmeans"].tolist())
    #print(df["total"].tolist())
    #print(df["kmeans"])
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #   print(df)
    if ksort == True:
        df = df.sort_values(by=["total", "kmeans"])

    if "kmeans" in  bar_lines:
        for i, x in enumerate(df["kmeans"].tolist()):
            if x > df["kmeans"].iloc[0]:
                ksplit = i
                break
    #if "length" in bar_lines:
    #   if length != None:
    #       for j, z in enumerate(df[col_name]):
    #           if z <= threshold:
    #               lsplit = j
    #               break
    #   else:
    #       lsplit = -10


    ## stacked bar plot ##
    if bar == True:
        xticks = np.arange(0, len(df))
        labels = df["sample"].tolist()#[i.replace("_filtered.vcf", "") for i in df["sample"].tolist()]
        fig, ax = plt.subplots()
        for i, v in enumerate(var):
            if i == 0:
                ax.bar(xticks, df[v].tolist(), 0.5, label = v, log=True)
            else:
                ax.bar(xticks, df[v].tolist(), 0.5, label = v, bottom=df[var[:i]].sum(axis=1), log=True)
        #df.plot(kind="bar", stacked=True)
        #plt.ylim([0, 1000])
        plt.tight_layout()
        plt.xlim([xticks[0]-1, xticks[-1]+1])
        plt.yscale(scale)
        plt.xticks(ticks = xticks, labels = labels, rotation = 90)
        plt.legend(var, loc = "upper left")
        if "kmeans" in  bar_lines:
            plt.axvline(x = ksplit-0.5, linestyle = "--", color="black")
        if "length" in bar_lines:
            plt.axvline(x = lsplit-0.5, linestyle = "--", color="red") 
        #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        #   print(df)
        
        if "kmeans" in  bar_numbers:
            for la, xp, yp in zip(df["kmeans"], xticks, df["total"]):
                plt.annotate(la, (xp-0.125, yp+1))
        
        if length != None:
            if bar_scale == True:
                print(max(df["total"]))
                sbh = max(df["total"])/20
                scale_bar = ax.barh(-5/2, max(xticks), sbh)
                gradientbars(scale_bar, df[col_name])
                norm=mpl.colors.Normalize(vmin=min(df[col_name]),vmax=max(df[col_name]))
                plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdYlGn"), ax=ax)
                plt.ylim(-5/2, max(df["total"])+3)
        if save == False:
            plt.show()
        else:
            plt.savefig(f"{file}_variant_type_stacked_bat_plot.png")
    
    ## box plot ##
    if box == True:
        if length != None:
            var_types = {"long": [], "short": []}
            mann_whitney = []
            is_bool = False
            if df[col_name].dtype == bool:
                is_bool = True
            var.append("total")
            for v in var:
                if is_bool == True:
                    var_types["long"].append(df[df[col_name] == False][v].tolist())
                    var_types["short"].append(df[df[col_name] == True][v].tolist())
                    #var_types[v] = [df[df[col_name] == False][v].tolist(), df[df[col_name] == True][v].tolist()]
                else:
                    var_types["long"].append(df[df["short"] == False][v].tolist())
                    var_types["short"].append(df[df["short"] == True][v].tolist())

            for i in range(len(var_types["long"])):
                mann_whitney.append(scipy.stats.mannwhitneyu(var_types["long"][i], var_types["short"][i]))
            #print(mann_whitney)
                    #var_types[v] = [df[df[col_name] < threshold][v].tolist(), df[df[col_name] >= threshold][v].tolist()]
            #print(var_types)
            lbox = plt.boxplot(var_types["long"], positions=np.array(range(len(var_types["long"])))*2.0-0.4, sym='', widths=0.6)
            plt.setp(lbox["boxes"], color = "blue")
            plt.setp(lbox["whiskers"], color = "blue")
            plt.setp(lbox["caps"], color = "blue")
            plt.setp(lbox["medians"], color = "blue")
            sbox = plt.boxplot(var_types["short"], positions=np.array(range(len(var_types["short"])))*2.0+0.4, sym='', widths=0.6)
            plt.setp(sbox["boxes"], color = "red")
            plt.setp(sbox["whiskers"], color = "red")
            plt.setp(sbox["caps"], color = "red")
            plt.setp(sbox["medians"], color = "red")
            plt.plot([], c="blue", label="long")
            plt.plot([], c="red", label="short")
            plt.legend()
            plt.xticks(ticks=range(0, len(var_types["short"])*2, 2), labels=var)
            for x, i in enumerate(mann_whitney):
                pval = i[1]
                if pval > 0.05:
                    continue
                elif 0.01 < pval <=0.5:
                    plt.annotate("*", (x*2, max(var_types["long"][x]+var_types["short"][x])), ha='center', fontsize="large")
                elif 0.001 < pval <= 0.01:
                    plt.annotate("**", (x*2, max(var_types["long"][x]+var_types["short"][x])), ha='center', fontsize="large")
                elif pval <= 0.001:
                    plt.annotate("***", (x*2, max(var_types["long"][x]+var_types["short"][x])), ha='center', fontsize="large")
            if save == False:
                plt.show()
            else:
                plt.savefig(f"{file}_variant_type_boxplot.png")
        else:
            var_types = {"lp": [], "p": []}
            mann_whitney = []
            var.append("total")
            for v in var:
                #var_types["lp"].append(get_len_array(df, v, "lowProb"))
                var_types["p"].append(get_len_array(df, v, "PASS"))
                #var_types[v] = [df[df[col_name] == False][v].tolist(), df[df[col_name] == True][v].tolist()]

            for i in range(len(var_types["lp"])):
                mann_whitney.append(scipy.stats.mannwhitneyu(var_types["lp"][i], var_types["p"][i]))
            #print(mann_whitney)
                    #var_types[v] = [df[df[col_name] < threshold][v].tolist(), df[df[col_name] >= threshold][v].tolist()]
            #print(var_types)
            lbox = plt.boxplot(var_types["lp"], positions=np.array(range(len(var_types["lp"])))*2.0-0.4, sym='', widths=0.6)
            plt.setp(lbox["boxes"], color = "blue")
            plt.setp(lbox["whiskers"], color = "blue")
            plt.setp(lbox["caps"], color = "blue")
            plt.setp(lbox["medians"], color = "blue")
            sbox = plt.boxplot(var_types["p"], positions=np.array(range(len(var_types["p"])))*2.0+0.4, sym='', widths=0.6)
            plt.setp(sbox["boxes"], color = "red")
            plt.setp(sbox["whiskers"], color = "red")
            plt.setp(sbox["caps"], color = "red")
            plt.setp(sbox["medians"], color = "red")
            plt.plot([], c="blue", label="lp")
            plt.plot([], c="red", label="p")
            plt.legend()
            plt.xticks(ticks=range(0, len(var_types["p"])*2, 2), labels=var)
            for x, i in enumerate(mann_whitney):
                pval = i[1]
                if pval > 0.05:
                    continue
                elif 0.01 < pval <=0.5:
                    plt.annotate("*", (x*2, max(var_types["lp"][x]+var_types["p"][x])), ha='center', fontsize="large")
                elif 0.001 < pval <= 0.01:
                    plt.annotate("**", (x*2, max(var_types["lp"][x]+var_types["p"][x])), ha='center', fontsize="large")
                elif pval <= 0.001:
                    plt.annotate("***", (x*2, max(var_types["lp"][x]+var_types["p"][x])), ha='center', fontsize="large")
            if save == False:
                plt.show()
            else:
                plt.savefig(f"{file}_variant_type_boxplot.png")


    ## histogram
    if hist == True:
        for v in var:
            plt.hist(df[v], alpha=0.5)
        plt.legend(var)
        if save == False:
            plt.show()
        else:
            plt.savefig(f"{file}_variant_type_histogram.png")


def gradientbars(bars, data):
    ax = bars[0].axes
    lim = ax.get_xlim()+ax.get_ylim()
    for bar in bars:      
        bar.set_zorder(1)
        bar.set_facecolor("none")
        x,y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        #grad = np.atleast_2d(np.linspace(min(data),1*w/max(data),256))
        ax.imshow(np.atleast_2d(data), extent=[x,x+w,y,y+h], aspect="auto", zorder=0, norm=mpl.colors.Normalize(vmin=min(data),vmax=max(data)), cmap="RdYlGn")
    ax.axis(lim)    
  


if __name__ == "__main__":
    plot_vars()
