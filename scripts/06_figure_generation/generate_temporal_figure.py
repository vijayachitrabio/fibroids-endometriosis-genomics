#!/usr/bin/env python3
"""
generate_temporal_figure.py
============================
Publication-quality butterfly dot plot for temporal architecture analysis.
Replaces the R-generated figure with proper trait labels and journal styling.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec

# ── Paths ────────────────────────────────────────────────────────────────────
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
DATA     = f"{repo_root}/results/Main_Table_2_Temporal_Full.csv"
OUT_DIR1 = f"{repo_root}/results"
OUT_DIR2 = f"{repo_root}/figures"
os.makedirs(OUT_DIR2, exist_ok=True)
FNAME    = "Figure_Temporal_Reviewer_Proof_Final.png"

# ── Trait label mapping (fix capitalisation) ─────────────────────────────────
TRAIT_LABELS = {
    "adenomyosis":           "Adenomyosis",
    "osteoarthritis":        "Osteoarthritis",
    "fibromyalgia":          "Fibromyalgia",
    "heavy_menstrual_bleeding": "Heavy Menstrual Bleeding",
    "endometrial_cancer":    "Endometrial Cancer",
    "ibs":                   "IBS",
    "nafld":                 "NAFLD",
    "endometrial_polyp":     "Endometrial Polyp",
    "depression":            "Depression",
    "anxiety":               "Anxiety",
    "endometriosis":         "Endometriosis",
    "osteoporosis":          "Osteoporosis",
    "iron_deficiency":       "Iron Deficiency Anaemia",
    "chronic_pelvic_pain":   "Chronic Pelvic Pain",
    "vitamin_d_deficiency":  "Vitamin D Deficiency",
    "ovarian_cancer":        "Ovarian Cancer",
    "cervical_cancer":       "Cervical Cancer",
    "pcos":                  "PCOS",
    "ptsd":                  "PTSD",
    "fibroids":              "Uterine Fibroids",
}

# ── Direction colours ─────────────────────────────────────────────────────────
DIR_COLOUR = {
    "Predominantly Post-Disease": "#1a3a6b",   # navy
    "Predominantly Pre-Disease":  "#c0392b",   # red
    "Mixed/Bidirectional":        "#7f8c8d",   # grey
}
DIR_SHORT = {
    "Predominantly Post-Disease": "Post-Disease",
    "Predominantly Pre-Disease":  "Pre-Disease",
    "Mixed/Bidirectional":        "Mixed",
}

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA)

# Exclude cancers as requested
cancers_to_exclude = ["endometrial_cancer", "ovarian_cancer", "cervical_cancer", "breast_cancer", "leiomyosarcoma"]
df = df[~df["trait"].isin(cancers_to_exclude)]

df["label"] = df["trait"].map(TRAIT_LABELS).fillna(df["trait"])

# ── Trait ordering ────────────────────────────────────────────────────────────
# Order by mean Prop_After across both diseases (ascending = pre-disease at bottom)
fib_df  = df[df["Disease"] == "Fibroids"].set_index("trait")
endo_df = df[df["Disease"] == "Endometriosis"].set_index("trait")

all_traits = sorted(
    set(fib_df.index) | set(endo_df.index),
    key=lambda t: (
        (fib_df.loc[t, "Prop_After_Conservative"]  if t in fib_df.index  else np.nan +
         endo_df.loc[t, "Prop_After_Conservative"] if t in endo_df.index else np.nan) / 2
        if (t in fib_df.index and t in endo_df.index)
        else (fib_df.loc[t, "Prop_After_Conservative"]  if t in fib_df.index  else
              endo_df.loc[t, "Prop_After_Conservative"] if t in endo_df.index else 0)
    )
)

# ── Sort traits by average Prop_After properly ───────────────────────────────
def avg_prop(t):
    vals = []
    if t in fib_df.index:
        vals.append(fib_df.loc[t, "Prop_After_Conservative"])
    if t in endo_df.index:
        vals.append(endo_df.loc[t, "Prop_After_Conservative"])
    return np.nanmean(vals) if vals else 0

all_traits = sorted(set(fib_df.index) | set(endo_df.index), key=avg_prop)
trait_labels = [TRAIT_LABELS.get(t, t) for t in all_traits]
n_traits = len(all_traits)
y_pos = {t: i for i, t in enumerate(all_traits)}

# ── Dot size scaling ──────────────────────────────────────────────────────────
def dot_size(n):
    """Scale N_total to marker area."""
    return max(55, min(500, n / 7))

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 12))
fig.patch.set_facecolor("white")

# Three columns: [endo panel] [labels] [fibroids panel]
gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 0.55, 1], wspace=0.04)
ax_endo = fig.add_subplot(gs[0, 0])   # Left: Endometriosis (x reversed)
ax_lab  = fig.add_subplot(gs[0, 1])   # Centre: trait labels
ax_fib  = fig.add_subplot(gs[0, 2])   # Right: Fibroids

# ── Helper: plot one panel ────────────────────────────────────────────────────
def plot_panel(ax, disease_df, reversed_x=False):
    for trait_key in all_traits:
        if trait_key not in disease_df.index:
            continue
        row = disease_df.loc[trait_key]
        yi   = y_pos[trait_key]
        prop = row["Prop_After_Conservative"]
        col  = DIR_COLOUR.get(row["Direction"], "#7f8c8d")
        sz   = dot_size(row["N_total"])
        high_sd = row["Flag_High_SameDay"]
        low_inf = row["Flag_Low_Inf"]

        if high_sd:
            ax.scatter(prop, yi, s=sz, c=col, marker="X",
                       zorder=4, linewidths=0.3, edgecolors="white")
        else:
            ax.scatter(prop, yi, s=sz, c=col, marker="o",
                       zorder=4, linewidths=0.3, edgecolors="white")
        # Flag low info with dashed outline
        if low_inf:
            ax.scatter(prop, yi, s=sz * 1.3, facecolors="none",
                       edgecolors=col, linewidths=1.2, linestyle="--", zorder=3)

# ── Plot endometriosis (x reversed) ──────────────────────────────────────────
plot_panel(ax_endo, endo_df, reversed_x=True)

ax_endo.set_xlim(0.98, -0.02)          # Reversed
ax_endo.set_ylim(-0.8, n_traits - 0.2)
ax_endo.set_xlabel("Proportion AFTER Diagnosis", fontsize=11, labelpad=8)
ax_endo.set_title("Endometriosis", fontsize=13, fontweight="bold", pad=10)
ax_endo.set_yticks([])
ax_endo.axvline(0.5, color="#999999", linestyle="--", linewidth=0.8, alpha=0.7)
ax_endo.set_xticks([0.0, 0.25, 0.50, 0.75])
ax_endo.set_xticklabels(["0.00", "0.25", "0.50", "0.75"], fontsize=9)
ax_endo.tick_params(axis="x", labelsize=9)
ax_endo.spines[["top", "right", "left"]].set_visible(False)
ax_endo.text(0.5, n_traits - 0.5, "← Post-disease", ha="center", va="top",
             fontsize=8, color="#777777", style="italic")

# Horizontal grid lines
for yi in range(n_traits):
    ax_endo.axhline(yi, color="#e8e8e8", linewidth=0.5, zorder=0)

# ── Plot fibroids ─────────────────────────────────────────────────────────────
plot_panel(ax_fib, fib_df, reversed_x=False)

ax_fib.set_xlim(-0.02, 0.98)
ax_fib.set_ylim(-0.8, n_traits - 0.2)
ax_fib.set_xlabel("Proportion AFTER Diagnosis", fontsize=11, labelpad=8)
ax_fib.set_title("Uterine Fibroids", fontsize=13, fontweight="bold", pad=10)
ax_fib.set_yticks([])
ax_fib.axvline(0.5, color="#999999", linestyle="--", linewidth=0.8, alpha=0.7)
ax_fib.set_xticks([0.0, 0.25, 0.50, 0.75])
ax_fib.set_xticklabels(["0.00", "0.25", "0.50", "0.75"], fontsize=9)
ax_fib.tick_params(axis="x", labelsize=9)
ax_fib.spines[["top", "left", "right"]].set_visible(False)
ax_fib.text(0.5, n_traits - 0.5, "Post-disease →", ha="center", va="top",
            fontsize=8, color="#777777", style="italic")

for yi in range(n_traits):
    ax_fib.axhline(yi, color="#e8e8e8", linewidth=0.5, zorder=0)

# ── Centre axis: trait labels ─────────────────────────────────────────────────
ax_lab.set_xlim(0, 1)
ax_lab.set_ylim(-0.8, n_traits - 0.2)
ax_lab.axis("off")

for trait_key, lbl in zip(all_traits, trait_labels):
    yi = y_pos[trait_key]
    ax_lab.text(0.5, yi, lbl, ha="center", va="center",
                fontsize=10, fontfamily="DejaVu Sans")

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = []
# Direction
for d, c in DIR_COLOUR.items():
    legend_elements.append(
        mlines.Line2D([0], [0], marker="o", color="none", markerfacecolor=c,
                      markeredgecolor="white", markersize=9,
                      label=DIR_SHORT[d])
    )
# Shape
legend_elements.append(
    mlines.Line2D([0], [0], marker="X", color="none", markerfacecolor="#555555",
                  markeredgecolor="white", markersize=9,
                  label="High same-day (>50%)")
)
# Dot sizes
for n_ex, lbl in [(500, "500"), (2000, "2,000"), (4000, "4,000")]:
    legend_elements.append(
        mlines.Line2D([0], [0], marker="o", color="none", markerfacecolor="#888888",
                      markeredgecolor="white",
                      markersize=np.sqrt(dot_size(n_ex)) * 0.75,
                      label=f"N = {lbl}")
    )

fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=4,
    fontsize=9,
    frameon=True,
    framealpha=0.95,
    edgecolor="#cccccc",
    bbox_to_anchor=(0.5, 0.005),
    handletextpad=0.5,
    columnspacing=1.2,
)

# ── Title ─────────────────────────────────────────────────────────────────────
fig.suptitle(
    "Temporal Architecture of Disease-Comorbidity Relationships\n"
    r"$\it{UK\ Biobank}$  ·  N ≥ 20 per cell  ·  Conservative same-day exclusion",
    fontsize=13, fontweight="bold", y=0.99, va="top"
)

fig.subplots_adjust(bottom=0.15, top=0.92)

# ── Save ──────────────────────────────────────────────────────────────────────
for out_dir in [OUT_DIR1, OUT_DIR2]:
    out_path = f"{out_dir}/{FNAME}"
    fig.savefig(out_path, dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_path}  ({os.path.getsize(out_path)//1000} KB)")

plt.close(fig)
print("Done.")
