#!/usr/bin/env python3
"""
plot_tissue_enrichment.py
=========================
Generate publication-quality 3-panel tissue enrichment figure from MAGMA
gene-property analysis results (.gsa.out files).

Input:   results/fibroids_tissue.gsa.out
         results/endometriosis_tissue.gsa.out
Output:  [Figures_Publication]/FigS_TissueEnrichment.png  (300 dpi)
         [Figures_Publication]/FigS_TissueEnrichment.tiff (600 dpi)
         [manuscript_prep_2026_04_17/figures]/FigS_TissueEnrichment.png

Usage:   python3 plot_tissue_enrichment.py
         (run from the magma_inputs/ directory, or edit BASE below)

Requirements:
    pandas, numpy, matplotlib  (all pre-installed in your environment)

Note: Run run_tissue_enrichment.sh first to generate the .gsa.out files.
"""

import os
import sys
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE       = os.path.dirname(SCRIPT_DIR)   # uterine_fibroids/
RESULTS    = os.path.join(SCRIPT_DIR, "results")
OUTDIR     = os.path.join(BASE, "outputs_v2", "Figures_Publication")
MANU       = os.path.join(BASE, "manuscript_prep_2026_04_17", "figures")
os.makedirs(OUTDIR, exist_ok=True)

# ── Configuration ─────────────────────────────────────────────────────────────
BONF_THRESH = 0.05          # Bonferroni significance (divided by n tissues)
FDR_THRESH  = 0.05          # FDR threshold (for colour coding)
TOP_N       = 54            # Show all 54 GTEx tissues (sorted by -log10P)

# Colour scheme — consistent with manuscript figure style
COL_SIG   = "#c0392b"   # red — Bonferroni significant
COL_FDR   = "#e67e22"   # orange — FDR significant but not Bonferroni
COL_NS    = "#aab7c4"   # grey — not significant
COL_REPRO = "#8e44ad"   # purple — reproductive / gynaecological tissues (highlight)

# Reproductive/gynaecological tissues to label regardless of significance
REPRO_TISSUES = {
    "Uterus", "Ovary", "Vagina", "Cervix_Ectocervix", "Cervix_Endocervix",
    "Fallopian_Tube", "Breast_Mammary_Tissue"
}

# Tissues to always display label if significant
ALWAYS_LABEL = {
    "Uterus", "Ovary", "Fallopian_Tube", "Vagina",
    "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes",
    "Whole_Blood", "Spleen", "Brain_Frontal_Cortex_BA9",
    "Brain_Cerebellum", "Artery_Aorta", "Heart_Left_Ventricle"
}


def load_gsa(trait):
    """Load and annotate a MAGMA .gsa.out file."""
    f = os.path.join(RESULTS, f"{trait}_tissue.gsa.out")
    if not os.path.exists(f):
        sys.exit(
            f"ERROR: {f} not found.\n"
            f"Please run run_tissue_enrichment.sh first."
        )
    df = pd.read_csv(f, sep=r'\s+', comment='#')
    # Expected columns: VARIABLE  TYPE  NGENES  BETA  BETA_STD  SE  P
    required = {"VARIABLE", "BETA", "P"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {f} missing columns: {missing}\nFound: {list(df.columns)}")

    n_tests         = len(df)
    df["neg_log10P"] = -np.log10(df["P"].clip(lower=1e-30))
    df["bonf"]       = df["P"] < (BONF_THRESH / n_tests)
    # FDR via Benjamini-Hochberg
    df_sorted        = df.sort_values("P").copy()
    ranks            = np.arange(1, len(df_sorted) + 1)
    df_sorted["fdr"] = (df_sorted["P"] * n_tests) / ranks
    # Make FDR monotone (step-up)
    fdr_vals = df_sorted["fdr"].values.copy()
    for i in range(len(fdr_vals) - 2, -1, -1):
        fdr_vals[i] = min(fdr_vals[i], fdr_vals[i + 1])
    df_sorted["fdr"] = fdr_vals
    df = df.merge(df_sorted[["VARIABLE", "fdr"]], on="VARIABLE", how="left")
    df["sig_level"] = "ns"
    df.loc[df["fdr"] < FDR_THRESH, "sig_level"] = "fdr"
    df.loc[df["bonf"], "sig_level"] = "bonf"
    df["is_repro"] = df["VARIABLE"].isin(REPRO_TISSUES)
    df["trait"]   = trait
    return df.sort_values("neg_log10P", ascending=False)


def draw_panel(ax, df, title, panel_label):
    """Draw a horizontal bar chart for one trait."""
    df_plot = df.sort_values("neg_log10P", ascending=True).copy()
    n       = len(df_plot)
    y_pos   = np.arange(n)

    # Bar colours
    colours = []
    for _, row in df_plot.iterrows():
        if row["is_repro"] and row["sig_level"] in ("bonf", "fdr"):
            colours.append(COL_REPRO)
        elif row["sig_level"] == "bonf":
            colours.append(COL_SIG)
        elif row["sig_level"] == "fdr":
            colours.append(COL_FDR)
        else:
            colours.append(COL_NS)

    ax.barh(y_pos, df_plot["neg_log10P"].values,
            color=colours, height=0.85, edgecolor="none", linewidth=0)

    # Significance threshold line
    n_tests = n
    bonf_line = -np.log10(BONF_THRESH / n_tests)
    ax.axvline(bonf_line, color="#2c3e50", linestyle="--", linewidth=0.9,
               label=f"Bonferroni (P={BONF_THRESH/n_tests:.3f})", zorder=5)

    # Y-axis tissue labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [t.replace("_", " ") for t in df_plot["VARIABLE"].values],
        fontsize=5.5
    )

    # Colour repro-tissue labels purple
    for tick, tissue in zip(ax.get_yticklabels(), df_plot["VARIABLE"].values):
        if tissue in REPRO_TISSUES:
            tick.set_color(COL_REPRO)
            tick.set_fontweight("bold")

    # X-axis
    max_val = max(df_plot["neg_log10P"].max() * 1.08, bonf_line * 1.3)
    ax.set_xlim(0, max_val)
    ax.set_xlabel(r"$-\log_{10}$ P value", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", length=0)

    # Top x-axis tick labels
    ax.xaxis.set_tick_params(labeltop=True, labelbottom=True)

    # Panel label (a, b, c)
    ax.text(-0.01, 1.02, panel_label, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="bottom", ha="right")

    # Title
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=6)

    # Clean spines
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)

    # Count label
    n_bonf = (df_plot["sig_level"] == "bonf").sum()
    n_fdr  = (df_plot["sig_level"] == "fdr").sum()
    ax.text(0.98, 0.02,
            f"Bonf: n={n_bonf}  FDR: n={n_fdr}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.5, color="#555555")


def make_shared_panel(fib, endo):
    """Compute shared-significance panel: tissues sig in BOTH traits."""
    merged = fib[["VARIABLE", "neg_log10P", "sig_level", "bonf", "is_repro"]].merge(
        endo[["VARIABLE", "neg_log10P", "sig_level", "bonf"]],
        on="VARIABLE", suffixes=("_fib", "_endo")
    )
    # Use geometric mean of -log10P for ranking
    merged["neg_log10P"]  = (merged["neg_log10P_fib"] + merged["neg_log10P_endo"]) / 2
    # sig_level = bonf only if bonf in BOTH
    merged["sig_level"]   = "ns"
    merged.loc[
        (merged["bonf_fib"]) | (merged["bonf_endo"]), "sig_level"] = "fdr"
    merged.loc[
        (merged["bonf_fib"]) & (merged["bonf_endo"]), "sig_level"] = "bonf"
    merged["bonf"]        = merged["bonf_fib"] & merged["bonf_endo"]
    merged["trait"]       = "shared"
    return merged.sort_values("neg_log10P", ascending=False)


# ── Load results ──────────────────────────────────────────────────────────────
print("Loading MAGMA gene-property results...")
fib  = load_gsa("fibroids")
endo = load_gsa("endometriosis")
shared = make_shared_panel(fib, endo)

print(f"  Fibroids:      {len(fib)} tissues  |  "
      f"Bonf sig: {fib['bonf'].sum()}  |  FDR sig: {(fib['sig_level']=='fdr').sum()}")
print(f"  Endometriosis: {len(endo)} tissues  |  "
      f"Bonf sig: {endo['bonf'].sum()}  |  FDR sig: {(endo['sig_level']=='fdr').sum()}")
print(f"  Shared:        Bonf in both: {shared['bonf'].sum()}")

# ── Figure layout ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(20, 18),
                          gridspec_kw={"wspace": 0.55})
fig.patch.set_facecolor("white")

draw_panel(axes[0], fib,    "Uterine fibroids\n(Gallagher 2019, N=218,728)",    "a")
draw_panel(axes[1], endo,   "Endometriosis\n(Rahmioglu 2023 / FinnGen R9)",     "b")
draw_panel(axes[2], shared, "Shared enrichment\n(significant in either trait)", "c")

# ── Shared legend ─────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=COL_SIG,   label="Bonferroni significant"),
    mpatches.Patch(color=COL_FDR,   label="FDR significant (BH)"),
    mpatches.Patch(color=COL_REPRO, label="Reproductive / gynaecological tissue"),
    mpatches.Patch(color=COL_NS,    label="Not significant"),
]
fig.legend(handles=legend_patches, fontsize=8,
           loc="lower center", bbox_to_anchor=(0.5, -0.01),
           ncol=4, frameon=True, framealpha=0.95, edgecolor="#cccccc")

# ── Supra-title ───────────────────────────────────────────────────────────────
fig.suptitle(
    "Tissue expression enrichment of GWAS signal — MAGMA gene-property analysis (GTEx v8)",
    fontsize=11, fontweight="bold", y=1.005
)
fig.text(0.5, 0.999,
         "54 GTEx tissues · gene-property model (direction=pos) · "
         "MAGMA v1.10 · 1000G EUR LD reference",
         ha="center", fontsize=7.5, color="#666666", style="italic")

# ── Save ──────────────────────────────────────────────────────────────────────
for fname, dpi in [("FigS_TissueEnrichment.png", 300),
                   ("FigS_TissueEnrichment.tiff", 600)]:
    out = os.path.join(OUTDIR, fname)
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}  ({os.path.getsize(out)//1000} KB)")

# Copy PNG to manuscript figures folder
manu_out = os.path.join(MANU, "FigS_TissueEnrichment.png")
shutil.copy(os.path.join(OUTDIR, "FigS_TissueEnrichment.png"), manu_out)
print(f"Copied to: {manu_out}")

plt.close(fig)
print("Done.")
