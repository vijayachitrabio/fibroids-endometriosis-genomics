"""Figure 4 (R1): cell-type heritability by S-LDSC, the two tests shown separately.

  A  heritability enrichment (baseline-LD v2.2) with SE; * = enrichment FDR < 0.05
  B  conditional coefficient Z (baseline v1.2);      * = coefficient FDR < 0.05

Former panels 4B/4C (curated pathway categories) were removed (Reviewer 1, comment 7).
FDR: Benjamini-Hochberg over the seven annotations within trait and model.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black", "axes.titlecolor": "black", "legend.labelcolor": "black"})
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

RR = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns")
OUT = RR.parent / "revision_figures"
C = {"Uterine fibroids": "#0072B2", "Endometriosis": "#D55E00"}
W = 0.38

d = pd.read_csv(RR / "SLDSC_full_statistics_from_raw_results.csv")
d["annotation"] = (d.annotation.str.replace("T_cell__CD8_", "CD8 T cell", regex=False)
                   .str.replace("Smooth_muscle", "Smooth muscle").str.replace("B_cell", "B cell"))
enr = d[d.baseline == "baseline-LD v2.2"]
cof = d[d.baseline == "baseline v1.2"]
order = (enr[enr.trait == "Uterine fibroids"].sort_values("enrichment").annotation.tolist())
y = np.arange(len(order))

fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6), dpi=600, gridspec_kw={"wspace": 0.12})
for ax, data, val, err, fdr, title, xlab, ref in (
        (axes[0], enr, "enrichment", "enrichment_SE", "enrichment_FDR",
         "A  Heritability enrichment (baseline-LD v2.2)", "Enrichment (proportion of h² / proportion of SNPs)", 1),
        (axes[1], cof, "coef_Z", None, "coef_FDR",
         "B  Conditional coefficient test (baseline v1.2)", "Coefficient Z score (τ / SE)", 0)):
    for k, (tr, sh) in enumerate((("Uterine fibroids", W / 2), ("Endometriosis", -W / 2))):
        t = data[data.trait == tr].set_index("annotation").loc[order]
        ax.barh(y + sh, t[val], W, color=C[tr], zorder=3, label=tr)
        if err:
            ax.errorbar(t[val], y + sh, xerr=t[err], fmt="none", ecolor="#333", elinewidth=0.7, capsize=2, zorder=4)
        for yi, (_, r) in zip(y, t.iterrows()):
            if r[fdr] < 0.05:
                x = r[val] + (r[err] if err else 0) + 0.12
                ax.text(max(x, 0.12), yi + sh, "*", va="center", fontsize=13, weight="bold", color="black")
    ax.axvline(ref, color="#4A5058", ls="--", lw=0.9, zorder=2)
    ax.set_title(title, loc="left", fontsize=11.5, weight="bold")
    ax.set_xlabel(xlab, fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#EDEFF2", lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=9.5)
axes[0].set_yticks(y)
axes[0].set_yticklabels(order)
axes[1].set_yticks(y)
axes[1].set_yticklabels([])
axes[0].set_xlim(0, None)
axes[1].text(0.99, 0.02, "* FDR < 0.05 (one-sided coefficient test)", transform=axes[1].transAxes,
             ha="right", fontsize=8.4, color="black")
axes[0].text(0.99, 0.02, "* FDR < 0.05 (enrichment test); bars ± SE", transform=axes[0].transAxes,
             ha="right", fontsize=8.4, color="black")
fig.legend(handles=[Line2D([0], [0], color=C[t], lw=7, label=t) for t in C], loc="lower center",
           ncol=2, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=0.1, right=0.985, top=0.92, bottom=0.17)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure_4_SLDSC_R1.{ext}", dpi=600, facecolor="white")
print("ok")
