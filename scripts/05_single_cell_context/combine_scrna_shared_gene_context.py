#!/usr/bin/env python3
"""Combine endometriosis and fibroid scRNA validation into one summary table."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "07_single_cell_context"
ASSETS = ROOT / "results" / "supplementary_tables"

ENDO = OUT / "Table_scEndoExplorer_Shared_MAGMA_CellType_Validation.csv"
FIBROID = OUT / "Table_GSE162122_Fibroid_scRNA_SharedGenes_CellContext.csv"

OUT_COMBINED = OUT / "Table_scRNA_Combined_Shared_MAGMA_KeyGene_Context.csv"
ASSET_COMBINED = ASSETS / "S14h_scRNA_Combined_Shared_MAGMA_KeyGene_Context.csv"


def top_fibroid_context(fibroid: pd.DataFrame, condition: str) -> pd.DataFrame:
    sub = fibroid[fibroid["Condition"] == condition].copy()
    major = sub[sub["Cells"] >= 1000].copy()
    if not major.empty:
        sub = major
    sub = sub.sort_values(["Gene", "Percent_expressing", "Mean_log1pCP10K"], ascending=[True, False, False])
    top = sub.groupby("Gene", as_index=False).first()
    return top.rename(
        columns={
            "Inferred_cell_type": f"{condition}_GSE162122_top_cell_context",
            "Mean_log1pCP10K": f"{condition}_GSE162122_mean_log1pCP10K",
            "Percent_expressing": f"{condition}_GSE162122_percent_expressing",
        }
    )[
        [
            "Gene",
            f"{condition}_GSE162122_top_cell_context",
            f"{condition}_GSE162122_mean_log1pCP10K",
            f"{condition}_GSE162122_percent_expressing",
        ]
    ]


def main() -> None:
    endo = pd.read_csv(ENDO)
    fibroid = pd.read_csv(FIBROID)

    endo_keep = endo[
        [
            "Gene",
            "MAGMA_log10P_fibroids",
            "MAGMA_log10P_endometriosis",
            "Endometriosis_scRNA_marker_cell_type",
            "Endometriosis_scRNA_marker_log2FC",
            "Endometriosis_scRNA_marker_FDR",
            "Endometriosis_celltype_LDSC_enrichment",
            "Endometriosis_celltype_LDSC_FDR",
            "Evidence_level",
        ]
    ].copy()

    combined = endo_keep.merge(top_fibroid_context(fibroid, "Fibroid"), on="Gene", how="left")
    combined = combined.merge(top_fibroid_context(fibroid, "Myometrium"), on="Gene", how="left")

    def interpretation(row):
        if row["Gene"] in {"ESR1", "GREB1"}:
            return "Shared hormone-responsive signal; endometriosis stromal support plus fibroid smooth-muscle/stromal expression."
        if row["Gene"] == "WNT4":
            return "Strong shared GWAS/MAGMA developmental locus; endometriosis stromal support, low but detectable fibroid expression."
        if row["Gene"] == "DNM3":
            return "Endothelial-context support in endometriosis and fibroid scRNA summaries."
        return "Supplementary context only; current scRNA support is weaker or nonspecific."

    combined["Manuscript_safe_interpretation"] = combined.apply(interpretation, axis=1)
    combined["Guardrail"] = (
        "GSE162122 fibroid cell contexts are marker-signature inferred from public count matrices; "
        "scEndoExplorer endometriosis contexts use precomputed markers/LDSC. Treat as validation/contextualization, not new primary discovery."
    )

    combined.to_csv(OUT_COMBINED, index=False)
    combined.to_csv(ASSET_COMBINED, index=False)
    print(f"Wrote {OUT_COMBINED}")
    print(f"Wrote {ASSET_COMBINED}")


if __name__ == "__main__":
    main()
