#!/usr/bin/env python3
"""
plot_magma_results.py
=====================
Run AFTER run_magma.sh completes.

Reads MAGMA gene-level results and produces:
  1. Figure_MAGMA_GeneLevel.png  — 4-panel figure replacing custom Fisher's analysis
  2. Table_MAGMA_TopGenes.csv    — Top 50 genes per trait for manuscript
  3. Supplementary_Table_S14a_fibroids_all_genes.tsv
  4. Supplementary_Table_S14b_endometriosis_all_genes.tsv
  5. Supplementary_Table_S14c_shared_genes.tsv
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

BASE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(BASE)
RESDIR = f"{BASE}/results"
FIGDIR = os.path.join(PROJECT_ROOT, "outputs_v2", "Figures_Publication")
SUPPDIR = os.path.join(PROJECT_ROOT, "outputs_v2")
GENE_LOC = os.path.join(BASE, "NCBI37.3.gene.loc")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(SUPPDIR, exist_ok=True)

BONFERRONI = 2.5e-6    # 0.05 / ~19,993 protein-coding genes
SUGGESTIVE = 1e-4

GENE_MAP = (
    pd.read_csv(
        GENE_LOC,
        sep="\t",
        header=None,
        names=["GENE", "CHR_LOC", "START_LOC", "STOP_LOC", "STRAND", "SYMBOL"],
    )[["GENE", "SYMBOL"]]
    .drop_duplicates()
)

# ── Load results ────────────────────────────────────────────────────────────
def load_magma(trait):
    f = f"{RESDIR}/{trait}_genes.genes.out"
    if not os.path.exists(f):
        sys.exit(f"ERROR: {f} not found. Run run_magma.sh first.")
    df = pd.read_csv(f, sep=r'\s+', comment="#")
    df = df.merge(GENE_MAP, on="GENE", how="left")
    df["SYMBOL"] = df["SYMBOL"].fillna(df["GENE"].astype(str))
    # Standard MAGMA columns: GENE SYMBOL NSNPS NPARAM N ZSTAT P
    df["LOG10P"] = -np.log10(df["P"].clip(lower=1e-300))
    df["TRAIT"]  = trait
    return df.sort_values("P")

print("Loading MAGMA results …")
fib  = load_magma("fibroids")
endo = load_magma("endometriosis")

# Merge for shared-gene analysis
merged = fib[["GENE","SYMBOL","LOG10P","P"]].merge(
    endo[["GENE","SYMBOL","LOG10P","P"]],
    on=["GENE","SYMBOL"], suffixes=("_fib","_endo")
)
merged["sig_fib"]  = merged["P_fib"]  < BONFERRONI
merged["sig_endo"] = merged["P_endo"] < BONFERRONI
merged["both_sig"] = merged["sig_fib"] & merged["sig_endo"]

shared     = merged[merged["both_sig"]].copy()
fib_only   = merged[merged["sig_fib"]  & ~merged["sig_endo"]].copy()
endo_only  = merged[merged["sig_endo"] & ~merged["sig_fib"]].copy()

print(f"  Fibroids Bonferroni-sig:      {merged['sig_fib'].sum():4d}")
print(f"  Endometriosis Bonferroni-sig: {merged['sig_endo'].sum():4d}")
print(f"  Shared:                       {len(shared):4d}")
print(f"  Fibroid-specific:             {len(fib_only):4d}")
print(f"  Endo-specific:                {len(endo_only):4d}")

# Colour classification for scatter
def classify(r):
    if r["both_sig"]:    return "shared"
    if r["sig_fib"]:     return "fib_only"
    if r["sig_endo"]:    return "endo_only"
    return "ns"

merged["class"] = merged.apply(classify, axis=1)
COLOURS = {"shared": "#1a3a6b", "fib_only": "#2a9d8f", "endo_only": "#e76f51", "ns": "#cccccc"}
SIZES   = {"shared": 30,        "fib_only": 20,        "endo_only": 20,        "ns": 4}
PANEL_DIR = os.path.join(FIGDIR, "Figure4_MAGMA_Panels")
os.makedirs(PANEL_DIR, exist_ok=True)

TITLE_FONT  = {"fontsize": 13, "fontweight": "bold"}
AXIS_FONT   = {"fontsize": 10}
TICK_FONT   = {"labelsize": 9}


def style_axis(ax):
    ax.tick_params(**TICK_FONT)
    ax.set_facecolor("#f9f9f9")


def maybe_set_title(ax, title):
    if title:
        ax.set_title(title, **TITLE_FONT)


def draw_panel_a(ax, title=None):
    top30_fib = fib.head(30).copy()
    top30_fib["colour"] = top30_fib["GENE"].isin(merged.loc[merged["both_sig"], "GENE"]).map(
        {True: "#1a3a6b", False: "#2a9d8f"}
    )
    ax.barh(
        range(len(top30_fib)),
        top30_fib["LOG10P"].values,
        color=top30_fib["colour"].values,
        edgecolor="white",
        linewidth=0.3,
    )
    ax.set_yticks(range(len(top30_fib)))
    ax.set_yticklabels(top30_fib["SYMBOL"].values, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(-np.log10(BONFERRONI), color="#cc0000", linestyle="--", linewidth=1, alpha=0.8)
    ax.set_xlabel("−log₁₀(P)", **AXIS_FONT)
    maybe_set_title(ax, title)
    style_axis(ax)
    leg_a = [Patch(facecolor="#1a3a6b", label="Shared signal"),
             Patch(facecolor="#2a9d8f", label="Fibroid-specific")]
    ax.legend(handles=leg_a, fontsize=8, loc="lower right")


def draw_panel_b(ax, title=None):
    top30_endo = endo.head(30).copy()
    top30_endo["colour"] = top30_endo["GENE"].isin(merged.loc[merged["both_sig"], "GENE"]).map(
        {True: "#1a3a6b", False: "#e76f51"}
    )
    ax.barh(
        range(len(top30_endo)),
        top30_endo["LOG10P"].values,
        color=top30_endo["colour"].values,
        edgecolor="white",
        linewidth=0.3,
    )
    ax.set_yticks(range(len(top30_endo)))
    ax.set_yticklabels(top30_endo["SYMBOL"].values, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(-np.log10(BONFERRONI), color="#cc0000", linestyle="--", linewidth=1, alpha=0.8)
    ax.set_xlabel("−log₁₀(P)", **AXIS_FONT)
    maybe_set_title(ax, title)
    style_axis(ax)
    leg_b = [Patch(facecolor="#1a3a6b", label="Shared signal"),
             Patch(facecolor="#e76f51", label="Endo-specific")]
    ax.legend(handles=leg_b, fontsize=8, loc="lower right")


def draw_panel_c(ax, title=None, annotate=True):
    for cls in ["ns", "fib_only", "endo_only", "shared"]:
        sub = merged[merged["class"] == cls]
        alpha = 0.35 if cls == "ns" else 0.75
        size = 3 if cls == "ns" else SIZES[cls]
        ax.scatter(
            sub["LOG10P_fib"],
            sub["LOG10P_endo"],
            c=COLOURS[cls],
            s=size,
            alpha=alpha,
            linewidths=0,
            label={"ns": "Non-significant", "fib_only": "Fibroid-specific",
                   "endo_only": "Endo-specific", "shared": "Shared"}[cls],
        )

    if annotate:
        shared_plot = merged[merged["both_sig"]].copy()
        shared_plot["sum_lp"] = shared_plot["LOG10P_fib"] + shared_plot["LOG10P_endo"]
        for _, row in shared_plot.nlargest(8, "sum_lp").iterrows():
            ax.annotate(
                row["SYMBOL"],
                (row["LOG10P_fib"], row["LOG10P_endo"]),
                fontsize=8,
                ha="left",
                va="bottom",
                xytext=(3, 2),
                textcoords="offset points",
            )

    ax.axvline(-np.log10(BONFERRONI), color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax.axhline(-np.log10(BONFERRONI), color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Fibroids −log₁₀(P)", **AXIS_FONT)
    ax.set_ylabel("Endometriosis −log₁₀(P)", **AXIS_FONT)
    maybe_set_title(ax, title)
    style_axis(ax)
    ax.legend(fontsize=8, markerscale=2, frameon=True, loc="upper left")


def draw_panel_d(ax, title=None):
    categories = [
        f"Shared\n(n={len(shared)})",
        f"Fibroid-specific\n(n={len(fib_only)})",
        f"Endo-specific\n(n={len(endo_only)})"
    ]
    counts  = [len(shared), len(fib_only), len(endo_only)]
    colours = ["#1a3a6b", "#2a9d8f", "#e76f51"]
    bars_d  = ax.bar(categories, counts, color=colours, edgecolor="white", width=0.55)
    for bar, cnt in zip(bars_d, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(cnt), ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_ylabel("Number of significant genes\n(Bonferroni P<2.5×10⁻⁶)", **AXIS_FONT)
    maybe_set_title(ax, title)
    style_axis(ax)
    ax.spines[["top","right"]].set_visible(False)


def save_single_panel(draw_fn, filename, figsize, **kwargs):
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    draw_fn(ax, **kwargs)
    fig.savefig(os.path.join(PANEL_DIR, filename), dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ── Figure ──────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor("white")
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

draw_panel_a(ax_a, title=None)
draw_panel_b(ax_b, title=None)
draw_panel_c(ax_c, annotate=True, title=None)
draw_panel_d(ax_d, title=None)

outfig = f"{FIGDIR}/Figure4_MAGMA_GeneLevel.png"
fig.savefig(outfig, dpi=180, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"\nFigure saved: {outfig}")

# Copy as supplementary too
import shutil
shutil.copy(outfig, f"{FIGDIR}/FigS7_MAGMA_Gene_Level.png")

save_single_panel(draw_panel_a, "Figure4A_MAGMA_Fibroids.png", (8, 9), title=None)
save_single_panel(draw_panel_b, "Figure4B_MAGMA_Endometriosis.png", (8, 9), title=None)
save_single_panel(draw_panel_c, "Figure4C_MAGMA_Scatter_Clean.png", (8, 7), annotate=False, title=None)
save_single_panel(draw_panel_d, "Figure4D_MAGMA_Overlap.png", (7, 6), title=None)
print(f"Individual panels: {PANEL_DIR}")

# ── Tables ───────────────────────────────────────────────────────────────────
# Main manuscript table: top 20 per trait side by side
top20_fib  = fib.head(20)[["SYMBOL","NSNPS","ZSTAT","P"]].reset_index(drop=True)
top20_endo = endo.head(20)[["SYMBOL","NSNPS","ZSTAT","P"]].reset_index(drop=True)
top20_fib.columns  = ["Gene_Fib","nSNPs_Fib","Z_Fib","P_Fib"]
top20_endo.columns = ["Gene_Endo","nSNPs_Endo","Z_Endo","P_Endo"]
table_main = pd.concat([top20_fib, top20_endo], axis=1)
table_main.to_csv(f"{SUPPDIR}/Table_MAGMA_TopGenes.csv", index=False)
print(f"Main table:  {SUPPDIR}/Table_MAGMA_TopGenes.csv")

# Supplementary tables
fib.to_csv(f"{SUPPDIR}/Supplementary_Table_S14a_Fibroids_MAGMA_AllGenes.tsv",
           sep="\t", index=False)
endo.to_csv(f"{SUPPDIR}/Supplementary_Table_S14b_Endometriosis_MAGMA_AllGenes.tsv",
            sep="\t", index=False)
shared.sort_values("P_fib").to_csv(
    f"{SUPPDIR}/Supplementary_Table_S14c_Shared_Genes_MAGMA.tsv",
    sep="\t", index=False)

print(f"Supp S14a:   {SUPPDIR}/Supplementary_Table_S14a_Fibroids_MAGMA_AllGenes.tsv")
print(f"Supp S14b:   {SUPPDIR}/Supplementary_Table_S14b_Endometriosis_MAGMA_AllGenes.tsv")
print(f"Supp S14c:   {SUPPDIR}/Supplementary_Table_S14c_Shared_Genes_MAGMA.tsv")
print("\nAll done.")
