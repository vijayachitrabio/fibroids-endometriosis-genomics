#!/usr/bin/env python3
"""Generate a temporal comorbidity network from the existing temporal table.

This script is additive: it reads existing CSV outputs and writes new figure
files only. It does not modify or delete any project data.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch


BASE = Path(__file__).resolve().parents[2]
INPUT = BASE / "results" / "Main_Table_2_Temporal_Full.csv"
PUB_DIR = BASE / "figures"
ASSET_DIR = BASE / "figures" / "main"

OUT_PNG = PUB_DIR / "Figure_Temporal_Network_Panel.png"
OUT_PDF = PUB_DIR / "Figure_Temporal_Network_Panel.pdf"
ASSET_PNG = ASSET_DIR / "Figure_3B_Temporal_Network_Panel.png"


LABELS = {
    "adenomyosis": "Adenomyosis",
    "anxiety": "Anxiety",
    "cervical_cancer": "Cervical cancer",
    "chronic_pelvic_pain": "Chronic pelvic pain",
    "depression": "Depression",
    "endometrial_cancer": "Endometrial cancer",
    "endometrial_polyp": "Endometrial polyp",
    "endometriosis": "Endometriosis",
    "fibroids": "Uterine fibroids",
    "fibromyalgia": "Fibromyalgia",
    "heavy_menstrual_bleeding": "Heavy menstrual bleeding",
    "ibs": "Irritable bowel syndrome",
    "iron_deficiency": "Iron deficiency",
    "nafld": "NAFLD",
    "osteoarthritis": "Osteoarthritis",
    "osteoporosis": "Osteoporosis",
    "ovarian_cancer": "Ovarian cancer",
    "pcos": "PCOS",
    "ptsd": "PTSD",
    "vitamin_d_deficiency": "Vitamin D deficiency",
}

COLORS = {
    "Predominantly Pre-Disease": "#f6c58f",
    "Predominantly Post-Disease": "#b9dce7",
    "Mixed/Bidirectional": "#d7cde8",
}


def pretty_trait(value: str) -> str:
    return LABELS.get(value, value.replace("_", " ").title())


def wrapped_label(value: str, width: int = 17) -> str:
    return "\n".join(textwrap.wrap(pretty_trait(value), width=width))


def edge_width(count: float, max_count: float) -> float:
    if max_count <= 0:
        return 0.8
    return 0.7 + 4.2 * math.sqrt(max(count, 0) / max_count)


def draw_arrow(ax, start, end, width, color="#111111", rad=0.0, alpha=0.96):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=16 + 2.2 * width,
        linewidth=width,
        color=color,
        alpha=alpha,
        shrinkA=18,
        shrinkB=22,
        connectionstyle=f"arc3,rad={rad}",
        zorder=2,
    )
    ax.add_patch(arrow)


def panel(ax, df: pd.DataFrame, disease: str, panel_letter: str):
    df = df.copy()
    df["trait_label"] = df["trait"].map(pretty_trait)
    df = df.sort_values(["Direction", "N_inf", "N_total"], ascending=[True, False, False])

    n = len(df)
    angles = np.linspace(math.radians(108), math.radians(468), n, endpoint=False)
    radius = 3.15
    positions = {
        row.trait: (radius * math.cos(angle), radius * math.sin(angle))
        for row, angle in zip(df.itertuples(index=False), angles)
    }
    center = (0.0, 0.0)
    max_node = max(float(df["N_total"].max()), 1.0)
    max_edge = max(float(df[["N_before", "N_after"]].to_numpy().max()), 1.0)

    # Edges first, so nodes and labels sit on top.
    for i, row in enumerate(df.itertuples(index=False)):
        pos = positions[row.trait]
        rad = 0.12 if i % 2 == 0 else -0.12
        if row.Direction == "Predominantly Pre-Disease":
            draw_arrow(ax, pos, center, edge_width(row.N_before, max_edge), rad=rad)
        elif row.Direction == "Predominantly Post-Disease":
            draw_arrow(ax, center, pos, edge_width(row.N_after, max_edge), rad=rad)
        else:
            draw_arrow(ax, pos, center, edge_width(row.N_before, max_edge), rad=0.16, alpha=0.72)
            draw_arrow(ax, center, pos, edge_width(row.N_after, max_edge), rad=-0.16, alpha=0.72)

    # Central disease node.
    central_count = int(df["N_inf"].sum())
    ax.scatter([0], [0], s=2300, color="#cf1212", edgecolors="white", linewidth=2.2, zorder=5)

    # Trait nodes and labels.
    for row in df.itertuples(index=False):
        x, y = positions[row.trait]
        size = 520 + 1550 * math.sqrt(float(row.N_total) / max_node)
        ax.scatter(
            [x],
            [y],
            s=size,
            color=COLORS.get(row.Direction, "#cccccc"),
            edgecolors="white",
            linewidth=1.5,
            zorder=4,
        )
        ha = "left" if x >= 0 else "right"
        dx = 0.18 if x >= 0 else -0.18
        ax.text(
            x + dx,
            y,
            wrapped_label(row.trait),
            ha=ha,
            va="center",
            fontsize=10.5,
            color="#111111",
            zorder=6,
        )

    ax.text(-3.75, -3.86, f"({panel_letter}) {disease}: temporal comorbidity network", fontsize=17, family="serif")
    ax.set_xlim(-4.05, 4.05)
    ax.set_ylim(-4.0, 4.0)
    ax.set_aspect("equal")
    ax.axis("off")


def main() -> None:
    df = pd.read_csv(INPUT)
    # Masked data handling for UKBB compliance
    for col in ["N_before", "N_after", "N_total", "N_inf", "N_same", "N_informative_non_same_day"]:
        if col in df.columns:
            df[col] = df[col].replace("<5", "1").astype(float)
            
    PUB_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=300)
    panel(axes[0], df[df["Disease"] == "Fibroids"], "Fibroids", "a")
    panel(axes[1], df[df["Disease"] == "Endometriosis"], "Endometriosis", "b")

    fig.suptitle(
        "Temporal clinical networks of comorbidity recording around gynecologic diagnosis",
        y=0.985,
        fontsize=18,
        weight="bold",
    )
    legend_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor=COLORS["Predominantly Pre-Disease"],
                   markeredgecolor="white", markersize=12, label="Mostly before index diagnosis"),
        plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor=COLORS["Predominantly Post-Disease"],
                   markeredgecolor="white", markersize=12, label="Mostly after index diagnosis"),
        plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor=COLORS["Mixed/Bidirectional"],
                   markeredgecolor="white", markersize=12, label="Mixed / bidirectional"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3, frameon=False, fontsize=11, bbox_to_anchor=(0.5, 0.02))
    fig.text(
        0.5,
        0.005,
        "Arrow direction indicates predominant temporal order; edge width scales with before/after count; node size scales with co-occurrence count.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.subplots_adjust(left=0.025, right=0.985, top=0.93, bottom=0.08, wspace=0.08)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(ASSET_PNG, dpi=300, bbox_inches="tight")
    print(OUT_PNG)
    print(OUT_PDF)
    print(ASSET_PNG)


if __name__ == "__main__":
    main()
