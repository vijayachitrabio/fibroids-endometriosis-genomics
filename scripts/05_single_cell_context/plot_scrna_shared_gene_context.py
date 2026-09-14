#!/usr/bin/env python3
"""Plot scRNA shared-gene cell-context summary for supplementary use."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "07_single_cell_context"
FIG_OUT = OUT / "Figures_Publication"
ASSET_FIG = ROOT / "figures" / "07_single_cell_context"

FIBROID = OUT / "Table_GSE162122_Fibroid_scRNA_SharedGenes_CellContext.csv"
ENDO = OUT / "Table_scEndoExplorer_Shared_MAGMA_CellType_Validation.csv"

GENES = ["ESR1", "GREB1", "WNT4", "DNM3"]
CELL_TYPES = ["Smooth_muscle", "Fibroblast_Stromal", "Endothelial", "Macrophage_Monocyte", "T_NK"]


def main() -> None:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    ASSET_FIG.mkdir(parents=True, exist_ok=True)

    fib = pd.read_csv(FIBROID)
    fib = fib[
        (fib["Condition"] == "Fibroid")
        & (fib["Gene"].isin(GENES))
        & (fib["Inferred_cell_type"].isin(CELL_TYPES))
    ].copy()
    fib["Gene"] = pd.Categorical(fib["Gene"], categories=GENES[::-1], ordered=True)
    fib["Inferred_cell_type"] = pd.Categorical(fib["Inferred_cell_type"], categories=CELL_TYPES, ordered=True)

    endo = pd.read_csv(ENDO)
    endo = endo[endo["Gene"].isin(GENES)].copy()

    fig, axes = plt.subplots(
        1, 2, figsize=(9.8, 4.6), dpi=600, gridspec_kw={"width_ratios": [1.35, 1]}
    )

    ax = axes[0]
    x = fib["Inferred_cell_type"].cat.codes
    y = fib["Gene"].cat.codes
    sizes = (fib["Percent_expressing"].clip(lower=0) + 0.2) * 22
    colours = fib["Mean_log1pCP10K"]
    sc = ax.scatter(x, y, s=sizes, c=colours, cmap="Greys", edgecolor="black", linewidth=0.35)
    ax.set_xticks(range(len(CELL_TYPES)))
    ax.set_xticklabels(
        ["Smooth\nmuscle", "Fibroblast/\nstromal", "Endothelial", "Macrophage/\nmonocyte", "T/NK"],
        fontsize=8,
        color="black",
    )
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES[::-1], fontsize=9, color="black", fontstyle="italic")
    ax.set_title("GSE162122 fibroid scRNA", fontsize=10, color="black")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(True, color="#e6e6e6", linewidth=0.5)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.6)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Mean log1p CP10K", fontsize=8, color="black")
    cbar.ax.tick_params(labelsize=7, colors="black")

    ax = axes[1]
    endo_plot = endo.dropna(subset=["Endometriosis_scRNA_marker_cell_type"]).copy()
    endo_plot["Gene"] = pd.Categorical(endo_plot["Gene"], categories=GENES[::-1], ordered=True)
    y2 = endo_plot["Gene"].cat.codes
    x2 = [0] * len(endo_plot)
    sizes2 = endo_plot["Endometriosis_scRNA_marker_log2FC"] * 115
    ax.scatter(x2, y2, s=sizes2, color="#666666", edgecolor="black", linewidth=0.35)
    for _, row in endo_plot.iterrows():
        ax.text(
            0.09,
            GENES[::-1].index(row["Gene"]),
            str(row["Endometriosis_scRNA_marker_cell_type"]),
            va="center",
            fontsize=8.5,
            color="black",
        )
    ax.set_xlim(-0.25, 1.15)
    ax.set_xticks([])
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES[::-1], fontsize=9, color="black", fontstyle="italic")
    ax.set_title("scEndoExplorer endometriosis markers", fontsize=10, color="black")
    ax.grid(True, axis="y", color="#e6e6e6", linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.6)

    fig.text(0.02, 0.96, "A", fontsize=11, weight="bold", color="black")
    fig.text(0.62, 0.96, "B", fontsize=11, weight="bold", color="black")
    fig.tight_layout(rect=[0.02, 0.03, 0.98, 0.95])

    for path in [
        FIG_OUT / "FigureS_scRNA_Shared_MAGMA_KeyGene_Context.png",
        ASSET_FIG / "FigureS_scRNA_Shared_MAGMA_KeyGene_Context.png",
    ]:
        fig.savefig(path, dpi=600, bbox_inches="tight")
    for path in [
        FIG_OUT / "FigureS_scRNA_Shared_MAGMA_KeyGene_Context.pdf",
        ASSET_FIG / "FigureS_scRNA_Shared_MAGMA_KeyGene_Context.pdf",
    ]:
        fig.savefig(path, bbox_inches="tight")
    print("Wrote scRNA context supplementary figure")


if __name__ == "__main__":
    main()
