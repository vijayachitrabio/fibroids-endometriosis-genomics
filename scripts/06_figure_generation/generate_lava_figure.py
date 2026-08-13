#!/usr/bin/env python3
"""
Generate publication-ready LAVA local genetic correlation figures.

Outputs:
  - Combined three-panel figure
  - Separate uncluttered panels A, B, and C as PNG and PDF

The script is intentionally self-contained and uses paths relative to the
current repository so it can be run from this workspace or from the GitHub copy.
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

DATA = REPO / "outputs_v2" / "LAVA_Local_rg_Results.csv"
FIGDIR = REPO / "outputs_v2" / "Figures_Publication"
MANDIR = REPO / "manuscript_prep_2026_04_17" / "figures"
GHDIR = REPO / "github_repo" / "figures"

for out_dir in (FIGDIR, MANDIR, GHDIR):
    out_dir.mkdir(parents=True, exist_ok=True)


COL_CONC = "#C43C30"
COL_DIV = "#1F4E79"
COL_NEUTRAL = "#7A7A7A"
COL_GRID = "#E6E6E6"
COL_BAND = "#F7F7F7"

CHROMS = list(range(1, 23)) + [23]
CHR_LABELS = [str(c) if c != 23 else "X" for c in CHROMS]

LABEL_GENES = {
    "WNT4": (0.95, 0.08),
    "STON2": (0.50, 0.11),
    "GREB1": (0.55, -0.14),
    "PDGFRA": (0.55, -0.12),
    "ESR1": (0.60, 0.10),
    "CDKN2B-AS1": (0.55, 0.10),
    "WT1": (0.55, 0.09),
    "HMGA2": (0.50, 0.10),
}


def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def load_lava():
    if not DATA.exists():
        raise FileNotFoundError(f"Could not find LAVA results: {DATA}")

    df = pd.read_csv(DATA)
    df = df[df["sig_FDR"] == True].copy()
    df["chrom"] = df["chrom"].astype(int)
    df["abs_rg"] = df["local_rg"].abs()
    df["direction"] = np.where(df["local_rg"] > 0, "Concordant", "Divergent")
    df["colour"] = np.where(df["local_rg"] > 0, COL_CONC, COL_DIV)

    chrom_max = df.groupby("chrom")["pos"].max().to_dict()

    def chrom_x(row):
        idx = CHROMS.index(row["chrom"])
        frac = row["pos"] / chrom_max.get(row["chrom"], 1)
        return idx + 0.10 + frac * 0.80

    df["xpos"] = df.apply(chrom_x, axis=1)
    df["dot_size"] = 22 + 150 * (df["abs_rg"] ** 2)

    fdr = pd.to_numeric(df["fdr_p"], errors="coerce")
    finite_pos = fdr[(fdr > 0) & np.isfinite(fdr)]
    floor = finite_pos.min() / 10 if len(finite_pos) else 1e-300
    df["fdr_plot"] = fdr.fillna(1.0).clip(lower=floor)
    df["neglog10_fdr"] = -np.log10(df["fdr_plot"])
    df["fdr_capped"] = fdr <= 0
    df["neglog10_fdr_cap"] = df.loc[~df["fdr_capped"], "neglog10_fdr"].max()
    if not np.isfinite(df["neglog10_fdr_cap"].iloc[0]):
        df["neglog10_fdr_cap"] = 300

    return df


def despine_grid(ax, axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, axis=axis, color=COL_GRID, linewidth=0.55, zorder=0)


def add_panel_label(ax, label, title):
    ax.set_title(f"{label}  {title}", loc="left", fontweight="bold", pad=8)


def gene_label_rows(df):
    rows = []
    for gene, (dx, dy) in LABEL_GENES.items():
        sub = df[df["nearest_gene"].fillna("") == gene]
        if not sub.empty:
            row = sub.loc[sub["abs_rg"].idxmax()].copy()
            row["label_gene"] = gene
            row["nudge_x"] = dx
            row["nudge_y"] = dy
            rows.append(row)

    chrx = df[(df["chrom"] == 23) & (df["local_rg"] < -0.40)]
    if not chrx.empty:
        row = chrx.loc[chrx["local_rg"].idxmin()].copy()
        row["label_gene"] = "chrX divergent"
        row["nudge_x"] = -1.00
        row["nudge_y"] = -0.05
        rows.append(row)

    return rows


def panel_a(ax, df, separate=False):
    n_conc = int((df["local_rg"] > 0).sum())
    n_div = int((df["local_rg"] < 0).sum())

    for i, _chrom in enumerate(CHROMS):
        if i % 2 == 0:
            ax.axvspan(i, i + 1, color=COL_BAND, zorder=0, linewidth=0)

    ax.axhline(0, color="#8C8C8C", linewidth=0.8, zorder=1)
    ax.axhline(0.5, color="#BDBDBD", linewidth=0.6, linestyle=":", zorder=1)

    for direction, color, zorder in (
        ("Divergent", COL_DIV, 2),
        ("Concordant", COL_CONC, 3),
    ):
        sub = df[df["direction"] == direction]
        ax.scatter(
            sub["xpos"], sub["local_rg"],
            c=color,
            s=sub["dot_size"] * (1.20 if separate else 1.0),
            alpha=0.82,
            linewidths=0.35,
            edgecolors="white",
            zorder=zorder,
        )

    for row in gene_label_rows(df):
        ax.annotate(
            row["label_gene"],
            xy=(row["xpos"], row["local_rg"]),
            xytext=(row["xpos"] + row["nudge_x"], row["local_rg"] + row["nudge_y"]),
            fontsize=9 if separate else 7.5,
            fontweight="bold",
            color=row["colour"],
            ha="center",
            va="bottom",
            arrowprops=dict(arrowstyle="-", color="#8C8C8C", lw=0.7),
        )

    ax.set_xticks([i + 0.5 for i in range(len(CHROMS))])
    ax.set_xticklabels(CHR_LABELS)
    ax.set_xlim(-0.30, len(CHROMS) + 0.10)
    ax.set_ylim(-0.78, 0.92)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Local genetic correlation (rg)")
    add_panel_label(
        ax,
        "A",
        f"Local genetic correlation across {len(df)} FDR-significant loci",
    )
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.1))
    despine_grid(ax, "y")

    handles = [
        mpatches.Patch(color=COL_CONC, label=f"Concordant rg > 0, n={n_conc}"),
        mpatches.Patch(color=COL_DIV, label=f"Divergent rg < 0, n={n_div}"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False)

    ax.text(
        0.995, 0.965,
        f"LAVA bivariate analysis\n2,495 genomic blocks\nFDR-significant loci: {len(df)}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5 if separate else 7.5,
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="#D0D0D0", alpha=0.94),
    )


def panel_b(ax, df, separate=False):
    n_conc = int((df["local_rg"] > 0).sum())
    n_div = int((df["local_rg"] < 0).sum())
    bins = np.linspace(-0.75, 0.85, 30)
    conc = df.loc[df["local_rg"] > 0, "local_rg"]
    div = df.loc[df["local_rg"] < 0, "local_rg"]

    ax.hist(
        conc, bins=bins, color=COL_CONC, alpha=0.82,
        edgecolor="white", linewidth=0.5, label=f"Concordant, n={n_conc}",
    )
    ax.hist(
        div, bins=bins, color=COL_DIV, alpha=0.82,
        edgecolor="white", linewidth=0.5, label=f"Divergent, n={n_div}",
    )

    ax.axvline(0, color="#4D4D4D", linewidth=1.0)
    median_conc = float(np.median(conc))
    median_all = float(np.median(df["local_rg"]))
    ax.axvline(median_conc, color=COL_CONC, linewidth=1.1, linestyle="--")
    ax.axvline(median_all, color=COL_NEUTRAL, linewidth=1.0, linestyle=":")

    ylim = ax.get_ylim()
    ax.text(
        median_conc + 0.025, ylim[1] * 0.92,
        f"Concordant median = {median_conc:.2f}",
        color=COL_CONC,
        fontsize=8.5 if separate else 7.3,
        va="top",
    )

    ax.set_xlabel("Local genetic correlation (rg)")
    ax.set_ylabel("Number of loci")
    add_panel_label(ax, "B", "Distribution of local rg")
    ax.legend(frameon=False, loc="upper left")
    despine_grid(ax, "y")


def panel_c(ax, df, separate=False):
    plot_df = df.copy()
    cap = float(plot_df["neglog10_fdr_cap"].iloc[0])
    ceiling = min(max(math.ceil(cap / 10) * 10, 40), 320)
    plot_df["evidence"] = plot_df["neglog10_fdr"].clip(upper=ceiling)
    plot_df["is_capped"] = plot_df["neglog10_fdr"] >= ceiling

    for direction, color, marker in (
        ("Divergent", COL_DIV, "o"),
        ("Concordant", COL_CONC, "o"),
    ):
        sub = plot_df[(plot_df["direction"] == direction) & (~plot_df["is_capped"])]
        ax.scatter(
            sub["local_rg"], sub["evidence"],
            c=color,
            s=sub["dot_size"] * (0.75 if separate else 0.55),
            alpha=0.76,
            linewidths=0.35,
            edgecolors="white",
            marker=marker,
            zorder=3,
            label=direction,
        )

        sub_cap = plot_df[(plot_df["direction"] == direction) & (plot_df["is_capped"])]
        ax.scatter(
            sub_cap["local_rg"], np.repeat(ceiling, len(sub_cap)),
            c=color,
            s=sub_cap["dot_size"] * (0.75 if separate else 0.55),
            alpha=0.82,
            linewidths=0.35,
            edgecolors="white",
            marker="^",
            zorder=4,
        )

    ax.axvline(0, color="#8C8C8C", linewidth=0.8)
    ax.axhline(-np.log10(0.05), color="#A0A0A0", linewidth=0.7, linestyle=":")

    label_candidates = []
    for gene in ("ESR1", "WT1", "WNT4", "GREB1", "CDKN2B-AS1"):
        sub = plot_df[plot_df["nearest_gene"].fillna("") == gene]
        if not sub.empty:
            label_candidates.append(sub.loc[sub["abs_rg"].idxmax()])
    chrx = plot_df[(plot_df["chrom"] == 23) & (plot_df["local_rg"] < -0.40)]
    if not chrx.empty:
        row = chrx.loc[chrx["local_rg"].idxmin()].copy()
        row["nearest_gene"] = "chrX divergent"
        label_candidates.append(row)

    offsets = {
        "ESR1": (-0.10, -22),
        "WT1": (-0.08, -48),
        "WNT4": (0.06, -30),
        "GREB1": (0.05, -18),
        "CDKN2B-AS1": (-0.04, 20),
        "chrX divergent": (0.04, -26),
    }

    for row in label_candidates:
        label = str(row["nearest_gene"])
        dx, dy = offsets.get(label, (0.04, -20))
        y = min(float(row["evidence"]), ceiling)
        ax.annotate(
            label,
            xy=(row["local_rg"], y),
            xytext=(row["local_rg"] + dx, max(4, y + dy)),
            fontsize=8.5 if separate else 7.0,
            fontweight="bold",
            color=row["colour"],
            ha="left" if dx >= 0 else "right",
            va="center",
            arrowprops=dict(arrowstyle="-", color="#8C8C8C", lw=0.65),
        )

    ax.set_xlim(-0.68, 0.84)
    ax.set_ylim(0, ceiling * 1.04)
    ax.set_xlabel("Local genetic correlation (rg)")
    ax.set_ylabel("-log10(FDR q-value)")
    add_panel_label(ax, "C", "Significance landscape")
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=COL_DIV,
                   label="Divergent", markersize=7),
        plt.Line2D([0], [0], marker="o", linestyle="", color=COL_CONC,
                   label="Concordant", markersize=7),
        plt.Line2D([0], [0], marker="^", linestyle="", color=COL_NEUTRAL,
                   label=f"Capped at -log10 FDR = {ceiling:.0f}", markersize=7),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    despine_grid(ax, "y")


def save_figure(fig, stem, dirs=(FIGDIR, MANDIR)):
    for out_dir in dirs:
        png = out_dir / f"{stem}.png"
        pdf = out_dir / f"{stem}.pdf"
        fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(pdf, bbox_inches="tight", facecolor="white")
        print(f"Saved: {png}")
        print(f"Saved: {pdf}")


def main():
    set_style()
    df = load_lava()
    n_conc = int((df["local_rg"] > 0).sum())
    n_div = int((df["local_rg"] < 0).sum())
    print(f"FDR-significant loci: {len(df)} (concordant={n_conc}, divergent={n_div})")

    fig = plt.figure(figsize=(15.5, 9.5))
    gs = GridSpec(
        2, 2, figure=fig, height_ratios=[1.45, 1.0],
        hspace=0.42, wspace=0.28, left=0.06, right=0.98, top=0.94, bottom=0.08,
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    panel_a(ax_a, df)
    panel_b(ax_b, df)
    panel_c(ax_c, df)
    save_figure(fig, "Figure_3_LAVA_Local_rg", dirs=(FIGDIR, MANDIR, GHDIR))
    plt.close(fig)

    panel_specs = [
        ("Figure_3A_LAVA_Local_rg_Chromosomal", (12.5, 4.9), panel_a),
        ("Figure_3B_LAVA_Local_rg_Distribution", (6.8, 4.8), panel_b),
        ("Figure_3C_LAVA_Local_rg_Significance", (7.2, 5.2), panel_c),
    ]
    for stem, size, func in panel_specs:
        fig, ax = plt.subplots(figsize=size)
        func(ax, df, separate=True)
        save_figure(fig, stem, dirs=(FIGDIR, MANDIR, GHDIR))
        plt.close(fig)

    print("Done.")


if __name__ == "__main__":
    main()
