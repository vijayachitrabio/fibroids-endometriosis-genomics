#!/usr/bin/env python3
"""
Paper 1, Figure 3: gene-level architecture.

  A  MAGMA fibroids       Bonferroni-significant genes ranked by -log10(P)
  B  MAGMA endometriosis  same

Only the ten genes significant in both disorders are labelled: labelling all
147 fibroid genes is unreadable, and the shared set is what the panel is for.
  C  cross-disorder comparison of gene-level significance

Counts reproduce the manuscript exactly: 147 fibroid genes, 38 endometriosis
genes, 10 significant in both, 137 significant only in the fibroid analysis,
28 only in the endometriosis analysis. Groups are threshold-defined (R1).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

BASE = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids")
TAB = BASE / "manuscript_assets_latest_2026_05_20" / "tables_supplementary"
OUT = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_figures")
IN_LAVA = {"SYNE1", "GREB1", "WNT4", "ESR1", "VIP", "LDLRAD2", "ARL14EP"}  # shared genes overlapping FDR-significant LAVA blocks (Table S2B)

THR = 2.5e-6                      # Bonferroni threshold used in the manuscript
C_SHARED = "#7B2D8E"              # significant in both
C_FIB = "#0072B2"                 # significant only in the fibroid analysis
C_ENDO = "#D55E00"                # significant only in the endometriosis analysis
N_LABEL = 22


def load():
    f = pd.read_csv(TAB / "S14a_MAGMA_Fibroids_AllGenes.tsv", sep="\t")
    e = pd.read_csv(TAB / "S14b_MAGMA_Endometriosis_AllGenes.tsv", sep="\t")
    for d in (f, e):
        d["SYMBOL"] = d["SYMBOL"].astype(str)
        d["nlp"] = -np.log10(d["P"].clip(lower=1e-320))
    return f, e


def panel_rank(ax, d, shared, colour, title, letter, n_top=6):
    sig = d[d["P"] < THR].sort_values("nlp", ascending=False).reset_index(drop=True)
    x = np.arange(len(sig))
    is_sh = sig["SYMBOL"].isin(shared).to_numpy()

    ax.vlines(x[~is_sh], 0, sig["nlp"][~is_sh], color=colour, lw=1.5, alpha=0.75)
    ax.vlines(x[is_sh], 0, sig["nlp"][is_sh], color=C_SHARED, lw=2.0)
    ax.scatter(x[~is_sh], sig["nlp"][~is_sh], s=13, color=colour, zorder=3)
    ax.scatter(x[is_sh], sig["nlp"][is_sh], s=30, color=C_SHARED, zorder=4,
               edgecolor="black", linewidth=0.4)

    # Rotated labels are taller than any sensible tier gap, so they are spread
    # horizontally on a single row instead, with slanted leaders back to the
    # true rank position.
    to_label = set(sig.head(n_top)["SYMBOL"]) | (set(sig["SYMBOL"]) & shared)
    lab = sig[sig["SYMBOL"].isin(to_label)].copy()
    lab["rank"] = [int(sig.index[sig["SYMBOL"] == g][0]) for g in lab["SYMBOL"]]
    lab = lab.sort_values("rank")

    ymax = sig["nlp"].max() * 1.62
    row = ymax * 0.70
    gap = len(sig) * 0.035
    cursor = -np.inf
    for _, r in lab.iterrows():
        xl = max(r["rank"], cursor + gap)
        cursor = xl
        sh = r["SYMBOL"] in shared
        ax.plot([r["rank"], xl], [r["nlp"] + ymax * 0.015, row - ymax * 0.012],
                color="#B9BFC7", lw=0.5, zorder=2)
        ax.annotate(r["SYMBOL"], (xl, row), rotation=90, fontsize=7.4,
                    ha="center", va="bottom",
                    color=C_SHARED if sh else "#22262A",
                    fontweight="bold" if sh else "normal", fontstyle="italic", zorder=6)

    ax.axhline(-np.log10(THR), color="#7A828A", ls="--", lw=0.8)
    ax.set_xlim(-2, len(sig) * 1.06)
    ax.set_ylim(0, sig["nlp"].max() * 1.62)
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=10)
    ax.set_xlabel(f"Bonferroni-significant genes, ranked (n = {len(sig)})", fontsize=9.5)
    ax.set_title(f"{letter}  {title}", loc="left", fontsize=11.5, weight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.grid(axis="y", color="#EDEFF2", lw=0.6)
    ax.set_axisbelow(True)


def panel_cross(ax, f, e, shared):
    m = f[["SYMBOL", "nlp"]].merge(e[["SYMBOL", "nlp"]], on="SYMBOL",
                                   suffixes=("_fib", "_endo"))
    thr = -np.log10(THR)
    sf = m["nlp_fib"] >= thr
    se = m["nlp_endo"] >= thr

    ax.scatter(m.loc[~sf & ~se, "nlp_fib"], m.loc[~sf & ~se, "nlp_endo"],
               s=5, color="#C9CFD6", alpha=0.55, linewidths=0, zorder=1)
    ax.scatter(m.loc[sf & ~se, "nlp_fib"], m.loc[sf & ~se, "nlp_endo"],
               s=20, color=C_FIB, alpha=0.85, linewidths=0, zorder=2)
    ax.scatter(m.loc[~sf & se, "nlp_fib"], m.loc[~sf & se, "nlp_endo"],
               s=20, color=C_ENDO, alpha=0.85, linewidths=0, zorder=2)
    both = m[sf & se]
    ax.scatter(both["nlp_fib"], both["nlp_endo"], s=52, color=C_SHARED,
               edgecolor="black", linewidth=0.5, zorder=4)
    inl = both[both.SYMBOL.isin(IN_LAVA)]
    ax.scatter(inl["nlp_fib"], inl["nlp_endo"], s=150, facecolor="none",
               edgecolor=C_SHARED, linewidth=1.1, zorder=3)

    # hand-placed labels (offset points) for the ten shared genes
    OFF = {"SYNE1": (8, 4), "GREB1": (8, 4), "ESR1": (8, 5), "ARL14EP": (-58, 8), "WNT4": (9, 5),
           "VIP": (-30, 9), "LDLRAD2": (-8, 15), "EFR3B": (10, 6), "DNM3": (12, -20), "PLEKHM3": (-62, -18)}
    for _, r in both.iterrows():
        dx, dy = OFF.get(r["SYMBOL"], (8, 4))
        ax.annotate(r["SYMBOL"], (r["nlp_fib"], r["nlp_endo"]),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=8.2, weight="bold", fontstyle="italic", color="#3A1145", zorder=5,
                    arrowprops=dict(arrowstyle="-", color="#9A8AA3", lw=0.5,
                                    shrinkA=0, shrinkB=3))

    ax.axvline(thr, color="#7A828A", ls="--", lw=0.8)
    ax.axhline(thr, color="#7A828A", ls="--", lw=0.8)
    ax.set_xlabel(r"Uterine fibroids  $-\log_{10}(P)$", fontsize=10)
    ax.set_ylabel(r"Endometriosis  $-\log_{10}(P)$", fontsize=10)
    ax.set_title("C  Cross-disorder comparison", loc="left",
                 fontsize=11.5, weight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.grid(color="#EDEFF2", lw=0.6)
    ax.set_axisbelow(True)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_SHARED,
               markeredgecolor="black", markersize=8,
               label=f"significant in both (n = {len(both)})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor=C_SHARED,
               markersize=12, label="shared gene in an FDR-significant LAVA block (n = 7)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_FIB,
               markersize=7, label=f"significant only in the fibroid analysis (n = {int((sf & ~se).sum())})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_ENDO,
               markersize=7, label=f"significant only in the endometriosis analysis (n = {int((~sf & se).sum())})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#C9CFD6",
               markersize=6, label="not significant"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.6)


def main():
    f, e = load()
    shared = set(f[f["P"] < THR]["SYMBOL"]) & set(e[e["P"] < THR]["SYMBOL"])
    print(f"fibroids {int((f['P'] < THR).sum())} · endometriosis "
          f"{int((e['P'] < THR).sum())} · shared {len(shared)}")

    fig = plt.figure(figsize=(13.6, 9.4), dpi=600)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.12], hspace=0.42, wspace=0.20)
    panel_rank(fig.add_subplot(gs[0, 0]), f, shared, C_FIB, "Uterine fibroids", "A", n_top=0)
    panel_rank(fig.add_subplot(gs[0, 1]), e, shared, C_ENDO, "Endometriosis", "B", n_top=0)
    panel_cross(fig.add_subplot(gs[1, :]), f, e, shared)

    fig.subplots_adjust(left=0.055, right=0.985, top=0.96, bottom=0.07)
    for ext in ("png", "pdf"):
        p = OUT / f"Figure_3_MAGMA_R1.{ext}"
        fig.savefig(p, dpi=600, facecolor="white")
        print("wrote", p)


if __name__ == "__main__":
    main()
