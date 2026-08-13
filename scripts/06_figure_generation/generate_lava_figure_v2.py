#!/usr/bin/env python3
"""
Generate LAVA figure matching the user's specific 3-panel style with error bars.
"""

from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parents[2]
if not (REPO / "outputs_v2").exists():
    REPO = Path.cwd()

DATA = REPO / "outputs_v2" / "formal_ld_audit" / "lava" / "Table_LAVA_Full_2495Blocks_AllResults_overlap_ldsc_corr_withFDR.csv"
FIGDIR = REPO / "Figures_600dpi_2026_08_03"
GHDIR = REPO / "github_repo" / "figures"

for out_dir in (FIGDIR, GHDIR):
    out_dir.mkdir(parents=True, exist_ok=True)

# Colours
COL_SIG = "#005AB5"
COL_NONSIG = "#808080"
COL_GRID = "#E0E0E0"

# Gene mappings
GENES = {
    "ESR1": (6, 151e6, 153e6),
    "WT1": (11, 29e6, 31e6),
    "WNT4": (1, 21e6, 23e6),
    "GREB1": (2, 10e6, 12e6),
}

def set_style():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

def load_data():
    df = pd.read_csv(DATA)
    df["chrom"] = df["CHR"].astype(int)
    df["sig"] = df["p_fdr"] < 0.05
    
    # Gene annotation
    df["gene"] = ""
    for gene, (chrom, start, end) in GENES.items():
        mask = (df["chrom"] == chrom) & (df["START"] <= end) & (df["STOP"] >= start)
        if mask.any():
            best_idx = df[mask]["p"].idxmin()
            df.at[best_idx, "gene"] = gene

    # Calculate positions for chromosome plot
    chrom_max = df.groupby("chrom")["STOP"].max().to_dict()
    chroms = list(range(1, 23))
    
    def chrom_x(row):
        c = row["chrom"]
        if c not in chroms: return -1
        idx = chroms.index(c)
        frac = row["START"] / chrom_max.get(c, 1)
        return idx + 1 + frac * 0.8 - 0.4  # Center around the integer tick
        
    df["xpos"] = df.apply(chrom_x, axis=1)
    df["neglog10_p"] = -np.log10(df["p"])
    
    return df, chroms

def plot_lollipops(ax, x, y, sig_mask):
    # Plot non-significant
    non_sig = ~sig_mask
    if non_sig.any():
        ax.vlines(x[non_sig], 0, y[non_sig], color=COL_NONSIG, linewidth=1.0, zorder=1)
        ax.plot(x[non_sig], y[non_sig], 'o', mfc='white', mec=COL_NONSIG, markersize=5, zorder=2)
        
    # Plot significant
    if sig_mask.any():
        ax.vlines(x[sig_mask], 0, y[sig_mask], color='#80B3E6', linewidth=1.5, zorder=3)
        ax.plot(x[sig_mask], y[sig_mask], 'o', mfc=COL_SIG, mec=COL_SIG, markersize=5, zorder=4)

def panel_a(ax, df, chroms):
    plot_lollipops(ax, df["xpos"], df["rho"], df["sig"])
    
    ax.axhline(0, color="#8C8C8C", linewidth=0.8, zorder=0)
    
    # Grid lines for chromosomes
    for c in chroms:
        ax.axvline(c, color=COL_GRID, linewidth=0.5, zorder=0)
        
    ax.set_xticks(chroms)
    ax.set_xticklabels([str(c) for c in chroms])
    ax.set_xlim(0.5, len(chroms) + 0.5)
    ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Local $r_g$")
    ax.set_title("A  Local genetic correlation by chromosome", loc="left", fontweight="bold")

def panel_b(ax, df):
    df_sorted = df.sort_values("rho").reset_index(drop=True)
    x = np.arange(1, len(df_sorted) + 1)
    
    plot_lollipops(ax, x, df_sorted["rho"], df_sorted["sig"])
    
    ax.axhline(0, color="#8C8C8C", linewidth=0.8, zorder=0)
    ax.axhline(1, color="#8C8C8C", linewidth=0.8, linestyle="--", zorder=0)
    
    ax.set_xlim(0, len(df_sorted) + 2)
    ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("Normally tested block (ranked)")
    ax.set_ylabel("Local $r_g$")
    ax.set_title("B  Ranked local genetic correlations", loc="left", fontweight="bold")

def panel_c(ax, df):
    sig = df["sig"]
    non_sig = ~sig
    
    ax.scatter(df.loc[non_sig, "rho"], df.loc[non_sig, "neglog10_p"], 
               facecolors='white', edgecolors=COL_NONSIG, s=20, zorder=1)
    ax.scatter(df.loc[sig, "rho"], df.loc[sig, "neglog10_p"], 
               facecolors=COL_SIG, edgecolors=COL_SIG, s=20, zorder=2)
               
    # Significance threshold line
    if sig.any():
        max_sig_p = df.loc[sig, "p"].max()
        thresh = -np.log10(max_sig_p)
        ax.axhline(thresh, color="#8C8C8C", linewidth=0.8, linestyle="--", zorder=0)
    
    # Annotate genes
    for _, row in df[df["gene"] != ""].iterrows():
        nudge_x = 0.03
        nudge_y = 2
        ha = "left"
        if row["rho"] > 0.8:
            nudge_x = -0.03
            ha = "right"
            
        ax.annotate(
            row["gene"],
            xy=(row["rho"], row["neglog10_p"]),
            xytext=(row["rho"] + nudge_x, row["neglog10_p"] + nudge_y),
            fontsize=7.5,
            fontweight="bold",
            color="black",
            ha=ha,
            va="bottom"
        )
        
    ax.set_xlim(-1.1, 1.1)
    # Add a bit of padding to top
    ax.set_ylim(-2, df["neglog10_p"].max() + 10)
    ax.set_xlabel("Local $r_g$")
    ax.set_ylabel("-$log_{10}(P)$")
    ax.set_title("C  Local correlation significance", loc="left", fontweight="bold")

def main():
    set_style()
    df, chroms = load_data()
    
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1], hspace=0.3, wspace=0.2)
    
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    
    panel_a(ax_a, df, chroms)
    panel_b(ax_b, df)
    panel_c(ax_c, df)
    
    out_png = FIGDIR / "Figure_3_LAVA_Local_rg_600dpi.png"
    out_pdf = FIGDIR / "Figure_3_LAVA_Local_rg_600dpi.pdf"
    
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    
    # Also save separate panels if needed? Not specifically requested, but the script only saves the main fig right now.
    print(f"Successfully generated LAVA figures at {out_png}")
    plt.close(fig)

if __name__ == "__main__":
    main()
