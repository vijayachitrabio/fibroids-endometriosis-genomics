#!/usr/bin/env python3
"""Supplementary Figure S1 (R1): single-cell context of shared MAGMA genes.

Copy of Figures_600dpi_2026_08_03/_support/generate_supplementary_s4.py with corrected labels only:
the dataset labelled "Endometriosis" is GSE218044 (adenomyosis uterine tissue, scEndoExplorer);
fibroid panels are GSE162122; the WT1 hub node is labelled by region. Data are unchanged."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
plt.rcParams.update({"text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black", "axes.titlecolor": "black", "legend.labelcolor": "black"})
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids")
SOURCE = ROOT / "source_data_2026_08_01/intermediate/source_data_panels.json"
OUTDIR = ROOT / "Funciton_int_ge/revision_figures"
PNG = OUTDIR / "Supplementary_Figure_S1_SingleCell_R1.png"
PDF = OUTDIR / "Supplementary_Figure_S1_SingleCell_R1.pdf"

TEXT = "black"
GRID = "#D9D9D9"
FIB = "#0072B2"
ENDO = "#C43C4A"
HORMONE = "#D33A2C"
SHARED = "#2474B5"
CONTEXT = "#2A9D68"
CORE = "#25364A"

GENES = ["ESR1", "GREB1", "WNT4", "DNM3"]
LABEL = {"Endometriosis": "Adenomyosis uterine tissue (GSE218044)", "Uterine fibroid": "Uterine fibroid tissue (GSE162122)"}
CELL_LABELS = {
    "Smooth_muscle": "Smooth\nmuscle",
    "Fibroblast_Stromal": "Fibroblast /\nstromal",
    "Endothelial": "Endothelial",
    "Macrophage_Monocyte": "Macrophage /\nmonocyte",
    "T_NK": "T / NK",
}


def load_sheet(payload: dict, name: str) -> pd.DataFrame:
    sheet = next(item for item in payload["sheets"] if item["name"] == name)
    return pd.DataFrame(sheet["rows"], columns=sheet["columns"])


def panel_label(ax, letter: str, x: float = -0.08) -> None:
    ax.text(x, 1.04, letter, transform=ax.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold", color=TEXT)


def style_axis(ax, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(TEXT)
    ax.spines["bottom"].set_color(TEXT)
    ax.tick_params(colors=TEXT, labelsize=7.0, width=0.7, length=3)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)


def plot_bubble(ax, data: pd.DataFrame) -> None:
    cell_types = list(CELL_LABELS)
    gene_to_y = {gene: len(GENES) - 1 - idx for idx, gene in enumerate(GENES)}
    cell_to_x = {cell: idx for idx, cell in enumerate(cell_types)}
    x = data["Inferred_cell_type"].map(cell_to_x).to_numpy()
    y = data["Gene"].map(gene_to_y).to_numpy()
    sizes = 11 + 38 * data["Percent_expressing"].to_numpy()
    scatter = ax.scatter(
        x,
        y,
        s=sizes,
        c=data["Mean_log1pCP10K"],
        cmap="Blues",
        vmin=0,
        vmax=max(0.20, float(data["Mean_log1pCP10K"].max())),
        edgecolors="#4A4A4A",
        linewidths=0.45,
    )
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels([CELL_LABELS[c] for c in cell_types], fontsize=6.5)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(list(reversed(GENES)), fontsize=7.3, fontstyle="italic")
    ax.set_xlim(-0.55, len(cell_types) - 0.45)
    ax.set_ylim(-0.55, len(GENES) - 0.45)
    ax.set_title("Uterine fibroid tissue (GSE162122)", loc="left", fontsize=9.0, fontweight="bold", color=TEXT, pad=6)
    ax.grid(color=GRID, linewidth=0.45)
    ax.set_axisbelow(True)

    size_handles = [
        ax.scatter([], [], s=11 + 38 * value, facecolors="white", edgecolors="#4A4A4A", linewidths=0.45, label=f"{value}%")
        for value in (1, 3, 6)
    ]
    legend = ax.legend(handles=size_handles, title="Expressing cells (%)", loc="lower center", bbox_to_anchor=(0.5, -0.22),
                       ncol=3, frameon=False, fontsize=6.3, title_fontsize=6.5, borderaxespad=0.0,
                       handletextpad=0.35, columnspacing=0.8)
    ax.add_artist(legend)
    cbar = ax.figure.colorbar(scatter, ax=ax, fraction=0.035, pad=0.035)
    cbar.set_label("Mean log1p CP10K", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6.1, width=0.5, length=2)
    panel_label(ax, "A")


def plot_scaled_context(ax, data: pd.DataFrame) -> None:
    x = np.arange(len(GENES))
    offsets = {"Endometriosis": -0.13, "Uterine fibroid": 0.13}
    styles = {
        "Endometriosis": dict(marker="D", facecolor=ENDO, edgecolor=ENDO),
        "Uterine fibroid": dict(marker="o", facecolor="white", edgecolor=FIB),
    }
    for dataset in ["Endometriosis", "Uterine fibroid"]:
        block = data[data["dataset"] == dataset].set_index("gene").loc[GENES]
        xx = x + offsets[dataset]
        for xpos, value in zip(xx, block["scaled_value"]):
            ax.vlines(xpos, 0, value, color=styles[dataset]["edgecolor"], linewidth=1.0)
        ax.scatter(xx, block["scaled_value"], s=45, linewidths=1.0, zorder=3, label=LABEL[dataset], **styles[dataset])
        for xpos, (_, row) in zip(xx, block.iterrows()):
            raw = f"{row.raw_value:.1f}%" if row.raw_value_unit == "percent expressing" else f"{row.raw_value:.1f}"
            align = "right" if dataset == "Endometriosis" else "left"
            ax.text(xpos, min(row.scaled_value + 0.065, 1.09), f"{row.dominant_context}\n{raw}", ha=align, va="bottom",
                    fontsize=5.7, color=TEXT, linespacing=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES, fontsize=7.3, fontstyle="italic")
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("Within-dataset scaled support", fontsize=7.5)
    ax.set_title("Dominant single-cell context", loc="left", fontsize=9.0, fontweight="bold", color=TEXT, pad=6)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=1, frameon=False, fontsize=6.4,
              handletextpad=0.4, columnspacing=1.1)
    style_axis(ax, "y")
    panel_label(ax, "B")


def plot_native_bar(ax, block: pd.DataFrame, dataset: str, colour: str, ylabel: str) -> None:
    block = block.set_index("gene").loc[GENES]
    x = np.arange(len(GENES))
    bars = ax.bar(x, block["raw_value"], color=colour, width=0.64)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES, fontsize=7.0, fontstyle="italic")
    ax.set_ylabel(ylabel, fontsize=7.2)
    ax.set_title(dataset, fontsize=8.2, fontweight="bold", color=TEXT, pad=5)
    upper = max(float(block["raw_value"].max()) * 1.30, 0.5)
    ax.set_ylim(0, upper)
    for bar, (_, row) in zip(bars, block.iterrows()):
        value = f"{row.raw_value:.1f}%" if row.raw_value_unit == "percent expressing" else f"{row.raw_value:.2f}"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + upper * 0.035,
                f"{row.dominant_context}\n{value}", ha="center", va="bottom", fontsize=5.6, color=TEXT, linespacing=1.0)
    style_axis(ax, "y")


def polar_xy(angle: float, radius: float) -> tuple[float, float]:
    radians = math.radians(angle)
    return radius * math.cos(radians), radius * math.sin(radians)


def plot_hub(ax, data: pd.DataFrame) -> None:
    nodes = data[data["record_type"] == "node"].copy()
    positions = {row.node_id: polar_xy(float(row.angle_degrees), float(row.radius)) for row in nodes.itertuples(index=False)}
    groups = {
        "Hormonal / reproductive hub": (0.95, 0.58),
        "Shared MAGMA genes": (-0.95, 0.02),
        "Fibroid structural context": (0.25, -0.78),
    }
    colours = {
        "Hormonal / reproductive hub": HORMONE,
        "Shared MAGMA genes": SHARED,
        "Fibroid structural context": CONTEXT,
    }

    for group, hub_pos in groups.items():
        linestyle = (0, (4, 3)) if group == "Fibroid structural context" else "solid"
        ax.plot([0, hub_pos[0]], [0, hub_pos[1]], color=colours[group], linewidth=1.0, linestyle=linestyle, alpha=0.55, zorder=1)
        for row in nodes[nodes["node_group"] == group].itertuples(index=False):
            x, y = positions[row.node_id]
            ax.plot([hub_pos[0], x], [hub_pos[1], y], color=colours[group], linewidth=0.9, linestyle=linestyle, alpha=0.48, zorder=1)

    ax.scatter([0], [0], s=600, color=CORE, edgecolors="white", linewidths=1.0, zorder=4)
    ax.text(0, 0, "Shared\ngenetic\ncore", ha="center", va="center", color="white", fontsize=6.0, fontweight="bold", linespacing=0.9, zorder=5)
    for group, (x, y) in groups.items():
        ax.scatter([x], [y], s=420, color=colours[group], edgecolors="white", linewidths=0.9, zorder=3)

    for row in nodes.itertuples(index=False):
        x, y = positions[row.node_id]
        colour = colours[row.node_group]
        is_shared = row.note == "shared MAGMA"
        ax.scatter([x], [y], s=float(row.marker_size) * 0.22, facecolors=colour if is_shared else "white",
                   edgecolors=colour, linewidths=1.2, zorder=3)
        norm = math.hypot(x, y)
        dx, dy = (x / norm * 0.24, y / norm * 0.24)
        if row.node_id == "STON2":
            dx, dy = 0.12, 0.13
        elif row.node_id == "WT1":
            dx, dy = 0.15, -0.10
        elif row.node_id == "HMGA2":
            dx, dy = 0.13, 0.10
        ha = "left" if x >= 0 else "right"
        if row.node_id == "WT1":
            ha, dx, dy = "center", 0.0, -0.36
        label = row.node_id if row.node_id != "WT1" else "chr11 30.2 Mb region"
        ax.text(x + dx, y + dy, label, ha=ha, va="center", fontsize=6.4, color=TEXT,
                fontstyle="italic" if row.node_id not in ("CDKN2B-AS1", "WT1") else "normal", zorder=5)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=HORMONE, markeredgecolor=HORMONE,
               markersize=5.5, label="Hormonal / reproductive shared MAGMA"),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=SHARED, markeredgecolor=SHARED,
               markersize=5.5, label="Other shared MAGMA"),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=CONTEXT,
               markersize=5.5, label="Fibroid-context gene"),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=HORMONE,
               markersize=5.5, label="LAVA-known locus"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.11), ncol=2, frameon=False,
              fontsize=5.7, handletextpad=0.3, columnspacing=0.8)
    ax.set_xlim(-3.85, 3.95)
    ax.set_ylim(-4.10, 3.20)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Conceptual shared and context-specific hubs", loc="left", fontsize=9.0, fontweight="bold", color=TEXT, pad=6)
    panel_label(ax, "D", x=-0.08)


def main() -> None:
    payload = json.loads(SOURCE.read_text())
    bubble = load_sheet(payload, "FigS4A")
    context = load_sheet(payload, "FigS4B")
    hub = load_sheet(payload, "FigS4D")
    if len(bubble) != 20 or len(context) != 8:
        raise ValueError("Unexpected S4 source-data dimensions")

    plt.rcParams.update({"font.family": "Arial", "text.color": TEXT, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.5, 8.2))
    outer = fig.add_gridspec(2, 2, width_ratios=[1.06, 0.94], height_ratios=[1.0, 1.02], hspace=0.38, wspace=0.34)
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    cgrid = outer[1, :].subgridspec(1, 2, wspace=0.25)
    ax_c1 = fig.add_subplot(cgrid[0, 0])
    ax_c2 = fig.add_subplot(cgrid[0, 1])

    plot_bubble(ax_a, bubble)
    plot_scaled_context(ax_b, context)
    plot_native_bar(ax_c1, context[context["dataset"] == "Endometriosis"], "Adenomyosis uterine tissue\n(GSE218044) markers", ENDO, "Marker log2FC")
    plot_native_bar(ax_c2, context[context["dataset"] == "Uterine fibroid"], "Uterine fibroid tissue\n(GSE162122)", FIB, "Cells expressing (%)")
    panel_label(ax_c1, "C", x=-0.12)
    panel_label(ax_c2, "D", x=-0.12)
    # Panel D (conceptual hub network) removed in R1: several nodes had no documented source.

    fig.subplots_adjust(left=0.07, right=0.985, top=0.965, bottom=0.08)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    metadata = {
        "Title": "Supplementary Figure S1. Single-cell context of selected genes significant in both MAGMA analyses",
        "Creator": "matplotlib",
    }
    fig.savefig(PNG, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.04, metadata=metadata)
    fig.savefig(PDF, facecolor="white", bbox_inches="tight", pad_inches=0.04, metadata=metadata)
    plt.close(fig)
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    main()
