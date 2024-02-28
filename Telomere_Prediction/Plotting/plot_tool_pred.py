import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from scipy.stats import *
from sklearn.metrics import mean_absolute_error

b = pd.read_csv("bgi_all_tool_pred.csv")
b = b.rename({"CompuTel": "computel"}, axis=1)
b["computel"] = b["computel"]
g = pd.read_csv("gel_cat_hunt_pred.csv")
g = g.rename({"telomerecat_pred": "telomerecat", "telomerehunter_pred": "telomerehunter"}, axis=1)

b["random"] = random.sample(range(int(min(b["stela"])*1000), int(max(b["stela"])*1000)), len(b))
b["random"] = b["random"]/1000
g["random"] = random.sample(range(int(min(g["stela"])*1000), int(max(g["stela"])*1000)), len(g))
g["random"] = g["random"]/1000

tools = ["telomerecat", "telomerehunter", "telseq", "computel", "qmotif", "random"]

bc = list(b.columns)
gc = list(g.columns)

fig, ax = plt.subplots(nrows=1, ncols=len(tools), constrained_layout=True)
fig.suptitle("Telomere prediction tools", fontsize="x-large")
for idx, t in enumerate(tools):
    bi = False
    gi = False
    if t in bc:
        btmp = b[["stela", t]]
        btmp = btmp.dropna()
        ax[idx].scatter(btmp["stela"], btmp[t], color="C0", alpha=0.5)
        bi = True
    if t in gc:
        gtmp = g[["stela", t]]
        gtmp = gtmp.dropna()
        ax[idx].scatter(gtmp["stela"], gtmp[t], color="C1", alpha=0.5) 
        gi = True
    lims = [
        np.min([ax[idx].get_xlim(), ax[idx].get_ylim()]),  # min of both ax[idx]es
        np.max([ax[idx].get_xlim(), ax[idx].get_ylim()]),  # max[idx] of both ax[idx]es
    ]
    if t == "telomerecat":
        lims = [0, 10]
    if t != "telseq":
        ax[idx].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax[idx].set_title(t)
    ax[idx].set_xlabel("stela (kb)")
    ax[idx].set_ylabel("prediction (kb)")
    ax[idx].set_aspect('equal')
    ax[idx].set_xlim(lims)
    ax[idx].set_ylim(lims)
    if bi and not gi:
        if t != "telseq":
            bmae = mean_absolute_error(btmp["stela"], btmp[t])
            bline = linregress(btmp["stela"], btmp[t])
            ax[idx].plot(b["stela"], bline.intercept+bline.slope*b["stela"], linestyle='--', color="C0")
            rleg = ax[idx].legend([f"{bline.rvalue:.3f}"], title="$R^2$", loc="upper right")
            rleg.legendHandles[-1].set_color("C0")
            ax[idx].add_artist(rleg)
            ax[idx].legend([f"{bmae:.3f}"], title="MAE:")
        if t == "telseq":
            gtmp = pd.read_csv("gel_telseq.csv")
            gtmp = gtmp.rename({"LENGTH_ESTIMATE": "telseq"}, axis=1)
            #gtmp = gtmp.drop_duplicates("stela").dropna()
            gtmp = gtmp.dropna()
            ax[idx].scatter(gtmp["stela"], gtmp[t], color="C1", alpha=0.5) 
            ax[idx].plot(lims, lims, 'k-', alpha=0.75, zorder=0)


            #bmae = mean_absolute_error(btmp["stela"], btmp[t])
            bline = linregress(btmp["stela"], btmp[t])
            gmae = mean_absolute_error(gtmp["stela"], gtmp[t])
            gline = linregress(gtmp["stela"], gtmp[t])
            a = pd.concat([btmp[["stela", t]], gtmp[["stela", t]]])
            a = a.dropna()
            aline = linregress(a["stela"], a[t])
            amae = mean_absolute_error(a["stela"], a[t])
            mleg = ax[idx].legend([f"{bmae:.3f}", f"{gmae:.3f}", f"{amae:.3f}"], title="MAE:", loc="upper left")
            mleg.legendHandles[0].set_color("C0")
            mleg.legendHandles[1].set_color("C1")
            mleg.legendHandles[2].set_color("black")
            ax[idx].add_artist(mleg)

            ax[idx].plot(b["stela"], bline.intercept+bline.slope*b["stela"], linestyle='--', color="C0")
            ax[idx].plot(g["stela"], gline.intercept+gline.slope*g["stela"], linestyle='--', color="C1")
            ax[idx].plot(a["stela"], aline.intercept+aline.slope*a["stela"], linestyle='--', color="C2")
            rleg = ax[idx].legend([f"{bline.rvalue:.3f}", f"{gline.rvalue:.3f}", f"{aline.rvalue:.3f}"], title="$R^2$", loc="upper right")
            rleg.legendHandles[0].set_color("C0")
            rleg.legendHandles[1].set_color("C1")
            rleg.legendHandles[2].set_color("C2")
            ax[idx].add_artist(rleg)
            ##ax[idx].legend([f"{bmae:.3f}", f"{gmae:.3f}", f"{amae:.3f}"], title="MAE:", loc="upper left")

    if gi and not bi:
        gmae = mean_absolute_error(gtmp["stela"], gtmp[t])
        gline = linregress(gtmp["stela"], gtmp[t])
        ax[idx].plot(g["stela"], gline.intercept+gline.slope*g["stela"], linestyle='--', color="C1")
        rleg = ax[idx].legend([f"{gline.rvalue:.3f}"], title="$R^2$", loc="upper right")
        rleg.legendHandles[-1].set_color("C1")
        ax[idx].add_artist(rleg)
        ax[idx].legend([f"{gmae:.3f}"], title="MAE:")
    if bi and gi:
        bmae = mean_absolute_error(btmp["stela"], btmp[t])
        bline = linregress(btmp["stela"], btmp[t])
        gmae = mean_absolute_error(gtmp["stela"], gtmp[t])
        gline = linregress(gtmp["stela"], gtmp[t])
        a = pd.concat([b[["stela", t]], g[["stela", t]]])
        a = a.dropna()
        aline = linregress(a["stela"], a[t])
        amae = mean_absolute_error(a["stela"], a[t])
        ax[idx].plot(b["stela"], bline.intercept+bline.slope*b["stela"], linestyle='--')
        ax[idx].plot(g["stela"], gline.intercept+gline.slope*g["stela"], linestyle='--')
        ax[idx].plot(a["stela"], aline.intercept+aline.slope*a["stela"], linestyle='--')
        rleg = ax[idx].legend([f"{bline.rvalue:.3f}", f"{gline.rvalue:.3f}", f"{aline.rvalue:.3f}"], title="$R^2$", loc="upper right")
        rleg.legendHandles[-1].set_color("C2")
        ax[idx].add_artist(rleg)
        ax[idx].legend([f"{bmae:.3f}", f"{gmae:.3f}", f"{amae:.3f}"], title="MAE:", loc="upper left")
fig.legend(["BGISeq 500", "Illumina HiSeq X/4000"],  loc="upper right", ncols=2, fontsize="large",  bbox_to_anchor=[0, 0, 1, 1]).set_zorder(10) 
plt.show()
