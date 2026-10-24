"""Figure 2 (R1): LAVA local genetic correlation with the R9 case/control counts.

  A  local rho by chromosome (88 valid tests), 95% CI, coloured by robustness
  B  the 16 FDR-significant blocks: rho with 95% CI under the primary
     sampling correlation (0.1724) and point estimates under 0.259 and 0.35
  C  rho versus -log10 P with regional (position) labels, not gene names

Addresses Reviewer 1 comment 2 and minor comment 3.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black", "axes.titlecolor": "black", "legend.labelcolor": "black"})
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

RR = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns")
OUT = RR.parent / "revision_figures"

C_CORE, C_SENS, C_NS = "#1F4E79", "#8DB3D9", "white"
STABLE = {"block_1074", "block_1075", "block_1076", "block_1077", "block_0024",
          "block_0216", "block_0666", "block_1642", "block_1644"}

loci = pd.read_csv(RR / "lava/lava_loci_2495_blocks.tsv", sep="\t")
v = pd.read_csv(RR / "lava/LAVA_corrected_N_ov01724_all_valid.csv")
sig = set(v[v.p_fdr < 0.05].LOC)
v["cls"] = np.where(v.LOC.isin(STABLE), "core", np.where(v.LOC.isin(sig), "sens", "ns"))
fill = {"core": C_CORE, "sens": C_SENS, "ns": C_NS}
size = {"core": 46, "sens": 46, "ns": 16}

chrlen = loci.groupby("CHR").STOP.max()
offset = chrlen.cumsum().shift(fill_value=0)
v["x"] = (v.START + v.STOP) / 2 + v.CHR.map(offset)


def region(r):
    return f"chr{r.CHR}: {r.START/1e6:.1f}–{r.STOP/1e6:.1f} Mb"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.grid(axis="y", color="#ECEEF1", lw=0.6)
    ax.set_axisbelow(True)


fig = plt.figure(figsize=(15.5, 10.6), dpi=600)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], width_ratios=[1.05, 1], hspace=0.27, wspace=0.42)

# ---------------------------------------------------------------- A
ax = fig.add_subplot(gs[0, :])
for c in ("ns", "sens", "core"):
    d = v[v.cls == c]
    ax.vlines(d.x, d["rho.lower"].clip(-1, 1), d["rho.upper"].clip(-1, 1),
              color="#C9CED4" if c == "ns" else fill[c], lw=0.5 if c == "ns" else 1.0, alpha=0.9, zorder=1)
    ax.scatter(d.x, d.rho, s=size[c], c=fill[c], edgecolor="black", linewidth=0.4, zorder=3)
for i, ch in enumerate(chrlen.index):
    if i % 2 == 0:
        ax.axvspan(offset[ch], offset[ch] + chrlen[ch], color="#F6F7F9", zorder=0)
ax.axhline(0, color="#7A828A", lw=0.7)
ax.set_xticks(offset + chrlen / 2)
ax.set_xticklabels(chrlen.index)
ax.set_xlim(0, (offset + chrlen).max())
ax.set_ylim(-1.08, 1.08)
ax.set_xlabel("Chromosome (GRCh37 LAVA blocks)", fontsize=10)
ax.set_ylabel("Local genetic correlation (ρ)", fontsize=10)
ax.set_title("A  Local genetic correlation across 88 blocks with valid bivariate tests", loc="left", fontsize=11.5, weight="bold")
style(ax)

# ---------------------------------------------------------------- B
ax = fig.add_subplot(gs[1, 0])
grid = {}
for k, f in {"0.259": "res_grid_corr_0.2590.csv", "0.35": "res_grid_corr_0.35.csv"}.items():
    d = pd.read_csv(RR / "lava" / f)
    grid[k] = d.set_index("LOC").rho
s16 = v[v.LOC.isin(sig)].copy()
s16["core"] = s16.LOC.isin(STABLE)
s16 = s16.sort_values(["core", "CHR", "START"], ascending=[False, True, True]).reset_index(drop=True)
y = np.arange(len(s16))[::-1]
for yi, (_, r) in zip(y, s16.iterrows()):
    col = C_CORE if r.core else C_SENS
    ax.hlines(yi, max(r["rho.lower"], -1), min(r["rho.upper"], 1), color=col, lw=2.2, zorder=2)
    ax.scatter(r.rho, yi, s=52, color=col, edgecolor="black", linewidth=0.5, zorder=4,
               marker=">" if r.rho >= 0.999 else "o")
    for k, mk in (("0.259", "s"), ("0.35", "D")):
        g = grid[k].get(r.LOC, np.nan)
        if pd.notna(g):
            ax.scatter(g, yi, s=22, marker=mk, facecolor="white", edgecolor="#333", linewidth=0.7, zorder=5)
ax.axvline(0, color="#7A828A", lw=0.7)
ax.axhline(y[s16.core.sum() - 1] - 0.5, color="#9AA1A9", lw=0.6, ls=":")
ax.text(-1.05, y[0] + 0.45, "Stable core", fontsize=9, weight="bold", color="black")
ax.text(-1.05, y[s16.core.sum()] + 0.25, "Sensitivity-dependent", fontsize=9, weight="bold", color="black")
ax.set_yticks(y)
ax.set_yticklabels([region(r) for _, r in s16.iterrows()], fontsize=8.6)
ax.set_xlim(-1.08, 1.08)
ax.set_ylim(-0.8, len(s16) + 0.6)
ax.set_xlabel("Local genetic correlation (ρ) with 95% CI", fontsize=10)
ax.set_title("B  FDR-significant blocks and sensitivity to sample overlap", loc="left", fontsize=11.5, weight="bold")
style(ax)
ax.grid(axis="y", visible=False)
ax.grid(axis="x", color="#ECEEF1", lw=0.6)
ax.legend(handles=[
    Line2D([0], [0], marker="o", color=C_CORE, mfc=C_CORE, mec="black", lw=2.2, ms=7, label="ρ (95% CI), sampling correlation 0.1724"),
    Line2D([0], [0], marker=">", color="none", mfc="#BBBBBB", mec="black", ms=7, label="estimate at the bound (ρ = 1)"),
    Line2D([0], [0], marker="s", color="none", mfc="white", mec="#333", ms=5.5, label="ρ, sampling correlation 0.259"),
    Line2D([0], [0], marker="D", color="none", mfc="white", mec="#333", ms=5, label="ρ, sampling correlation 0.35"),
], loc="lower left", fontsize=8, frameon=True, framealpha=0.95, edgecolor="#DDD")

# ---------------------------------------------------------------- C
ax = fig.add_subplot(gs[1, 1])
v["nlp"] = -np.log10(v.p.clip(lower=1e-300))
for c in ("ns", "sens", "core"):
    d = v[v.cls == c]
    ax.scatter(d.rho, d.nlp, s=size[c] + 6, c=fill[c], edgecolor="black", linewidth=0.4, zorder=3)
fdr_line = v[v.p_fdr < 0.05].p.max()
ax.axhline(-np.log10(fdr_line), color="#7A828A", ls="--", lw=0.8)
ax.text(-1.0, -np.log10(fdr_line) + 1.2, "FDR = 0.05", fontsize=8, color="black")
labs = {"block_1075": "chr6: 152.4–153.1 Mb", "block_0024": "chr1: 22.1–22.7 Mb", "block_1076": "chr6: 153.1–153.8 Mb",
        "block_1074": "chr6: 151.6–152.4 Mb", "block_0216": "chr2: 11.2–12.0 Mb", "block_1642": "chr11: 29.5–30.6 Mb"}
offs = {"block_1075": (-115, -4), "block_0024": (-100, 10), "block_1076": (-125, 20), "block_1074": (-150, 3),
        "block_0216": (-165, -6), "block_1642": (-150, -7)}
for loc, lab in labs.items():
    r = v[v.LOC == loc].iloc[0]
    ax.annotate(lab, (r.rho, r.nlp), xytext=offs[loc], textcoords="offset points", fontsize=8.4,
                arrowprops=dict(arrowstyle="-", color="#888", lw=0.5))
ax.set_xlim(-1.08, 1.08)
ax.set_xlabel("Local genetic correlation (ρ)", fontsize=10)
ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=10)
ax.set_title("C  Strength of local association", loc="left", fontsize=11.5, weight="bold")
style(ax)

fig.legend(handles=[
    Line2D([0], [0], marker="o", color="none", mfc=C_CORE, mec="black", ms=8, label="FDR < 0.05, stable across sampling correlations (9 blocks, 5 regions)"),
    Line2D([0], [0], marker="o", color="none", mfc=C_SENS, mec="black", ms=8, label="FDR < 0.05, sensitivity-dependent (7 blocks)"),
    Line2D([0], [0], marker="o", color="none", mfc="white", mec="black", ms=5, label="not significant"),
], loc="lower center", ncol=3, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
fig.subplots_adjust(left=0.1, right=0.985, top=0.965, bottom=0.085)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure_2_LAVA_R1.{ext}", dpi=600, facecolor="white")
print("ok", len(v), len(sig), fdr_line)
