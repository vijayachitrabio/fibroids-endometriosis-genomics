#!/usr/bin/env python3
"""Add exact temporal-ordering tests to the existing temporal table.

Reads the current temporal summary and writes new supplementary outputs only.
Existing data, figures, and scripts are not modified.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest


BASE = Path(__file__).resolve().parents[2]
INPUT = BASE / "results" / "Main_Table_2_Temporal_Full.csv"
OUT = BASE / "results" / "Table_Temporal_Statistical_Tests.csv"
ASSET_OUT = BASE / "results" / "S11b_Temporal_Statistical_Tests.csv"


def bh_fdr(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = ~np.isnan(p)
    vals = p[mask]
    m = len(vals)
    if m == 0:
        return pd.Series(out, index=p_values.index)
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty(m, dtype=float)
    tmp[order] = adj
    out[mask] = tmp
    return pd.Series(out, index=p_values.index)


def exact_p(after: int, before: int) -> float:
    n = after + before
    if n <= 0:
        return np.nan
    return float(binomtest(after, n=n, p=0.5, alternative="two-sided").pvalue)


def exact_direction(after: int, before: int, fdr: float) -> str:
    if np.isnan(fdr) or fdr >= 0.05:
        return "Not significant"
    if after > before:
        return "Post-diagnosis enriched"
    if before > after:
        return "Pre-diagnosis enriched"
    return "Balanced"


def main() -> None:
    df = pd.read_csv(INPUT)

    # Primary test: ignore same-day administrative co-coding and compare
    # post- vs pre-index dated events.
    df["N_informative_non_same_day"] = df["N_before"] + df["N_after"]
    df["Prop_After_NonSameDay"] = df["N_after"] / df["N_informative_non_same_day"]
    df["Temporal_Binomial_P"] = [
        exact_p(int(a), int(b)) for a, b in zip(df["N_after"], df["N_before"])
    ]
    df["Temporal_Binomial_FDR"] = df.groupby("Disease", group_keys=False)["Temporal_Binomial_P"].apply(bh_fdr)
    df["Temporal_Test_Direction"] = [
        exact_direction(int(a), int(b), float(q))
        for a, b, q in zip(df["N_after"], df["N_before"], df["Temporal_Binomial_FDR"])
    ]

    # Sensitivity 1: conservative post-diagnosis proportion already used in
    # the manuscript, retaining same-day events in the denominator.
    df["Sensitivity_Conservative_Post_Prop"] = df["N_after"] / df["N_total"]
    df["Sensitivity_Conservative_Class"] = np.select(
        [
            df["Sensitivity_Conservative_Post_Prop"] > 0.60,
            df["Sensitivity_Conservative_Post_Prop"] < 0.40,
        ],
        ["Post-dominant", "Pre/same-day dominant"],
        default="Mixed",
    )

    # Sensitivity 2: stress-test ambiguity by assigning same-day events to
    # either side. If both extremes remain on the same side of 0.5, the
    # temporal direction is robust to same-day handling.
    df["Sensitivity_Post_Prop_SameDay_As_Pre"] = df["N_after"] / df["N_total"]
    df["Sensitivity_Post_Prop_SameDay_As_Post"] = (df["N_after"] + df["N_same"]) / df["N_total"]
    robust_post = df["Sensitivity_Post_Prop_SameDay_As_Pre"] > 0.5
    robust_pre = df["Sensitivity_Post_Prop_SameDay_As_Post"] < 0.5
    df["Sensitivity_SameDay_Robust_Direction"] = np.select(
        [robust_post, robust_pre],
        ["Robust post-diagnosis", "Robust pre-diagnosis"],
        default="Same-day-sensitive / ambiguous",
    )

    df["High_SameDay_Fraction"] = df["N_same"] / df["N_total"]

    ordered = [
        "Disease", "trait", "N_total", "N_before", "N_after", "N_same",
        "N_informative_non_same_day", "Prop_After_NonSameDay",
        "Temporal_Binomial_P", "Temporal_Binomial_FDR", "Temporal_Test_Direction",
        "Sensitivity_Conservative_Post_Prop", "Sensitivity_Conservative_Class",
        "Sensitivity_Post_Prop_SameDay_As_Pre", "Sensitivity_Post_Prop_SameDay_As_Post",
        "Sensitivity_SameDay_Robust_Direction", "High_SameDay_Fraction",
        "Med_Lag_Post_Yrs", "Med_Lag_Pre_Yrs", "Flag_High_SameDay", "Direction",
    ]
    out = df[ordered].sort_values(["Disease", "Temporal_Binomial_FDR", "trait"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    ASSET_OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    out.to_csv(ASSET_OUT, index=False)

    summary = (
        out.assign(sig=out["Temporal_Binomial_FDR"] < 0.05)
        .groupby("Disease")
        .agg(
            traits=("trait", "count"),
            fdr_significant=("sig", "sum"),
            post_enriched=("Temporal_Test_Direction", lambda s: int((s == "Post-diagnosis enriched").sum())),
            pre_enriched=("Temporal_Test_Direction", lambda s: int((s == "Pre-diagnosis enriched").sum())),
            same_day_sensitive=("Sensitivity_SameDay_Robust_Direction", lambda s: int((s == "Same-day-sensitive / ambiguous").sum())),
        )
    )
    print(OUT)
    print(ASSET_OUT)
    print(summary.to_string())


if __name__ == "__main__":
    main()
