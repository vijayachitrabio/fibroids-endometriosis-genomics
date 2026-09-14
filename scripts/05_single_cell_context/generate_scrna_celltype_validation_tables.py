#!/usr/bin/env python3
"""Create additive scRNA cell-type validation tables for shared MAGMA genes.

Current inputs:
- scEndoExplorer endometriosis marker genes and LDSC cell-type enrichment.
- Local MAGMA shared-gene table from this project.

The fibroid GSE162122 parser can be added after the processed GEO matrix is
available; this script intentionally does not modify any existing inputs.
"""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "07_single_cell_context"
OUT = ROOT / "results" / "07_single_cell_context"
ASSETS = ROOT / "results" / "supplementary_tables"

SHARED_GENES = DATA / "Supplementary_Table_S14c_Shared_Genes_MAGMA.tsv"
ENDO_MARKERS = DATA / "cell_type_marker_genes.csv"
ENDO_LDSC = DATA / "ldsc_cell_type_enrichment.csv"

OUT_ENDO = OUT / "Table_scEndoExplorer_Shared_MAGMA_CellType_Validation.csv"
ASSET_ENDO = ASSETS / "S14f_scEndoExplorer_Shared_MAGMA_CellType_Validation.csv"


def load_shared_genes() -> pd.DataFrame:
    shared = pd.read_csv(SHARED_GENES, sep="\t")
    return shared.rename(
        columns={
            "SYMBOL": "Gene",
            "LOG10P_fib": "MAGMA_log10P_fibroids",
            "LOG10P_endo": "MAGMA_log10P_endometriosis",
            "P_fib": "MAGMA_P_fibroids",
            "P_endo": "MAGMA_P_endometriosis",
        }
    )


def endometriosis_marker_summary(shared: pd.DataFrame) -> pd.DataFrame:
    markers = pd.read_csv(ENDO_MARKERS)
    ldsc = pd.read_csv(ENDO_LDSC)

    marker_hits = markers[markers["gene"].isin(shared["Gene"])].copy()
    marker_hits = marker_hits.sort_values(["gene", "p_val_adj", "avg_log2FC"])

    records = []
    for row in shared.itertuples(index=False):
        hits = marker_hits[marker_hits["gene"] == row.Gene]
        if hits.empty:
            records.append(
                {
                    "Gene": row.Gene,
                    "Endometriosis_scRNA_marker_cell_type": "",
                    "Endometriosis_scRNA_marker_log2FC": None,
                    "Endometriosis_scRNA_marker_FDR": None,
                    "Endometriosis_all_marker_hits": "",
                    "Endometriosis_celltype_LDSC_enrichment": "",
                    "Endometriosis_celltype_LDSC_FDR": None,
                    "Evidence_level": "No marker hit in scEndoExplorer marker table",
                }
            )
            continue

        best = hits.iloc[0]
        ldsc_hit = ldsc[ldsc["cell_type"] == best["cell_type"]]
        ldsc_enrichment = None
        ldsc_fdr = None
        if not ldsc_hit.empty:
            ldsc_enrichment = float(ldsc_hit.iloc[0]["enrichment"])
            ldsc_fdr = float(ldsc_hit.iloc[0]["fdr"])

        all_hits = "; ".join(
            f"{h.cell_type} log2FC={h.avg_log2FC:.2f}, FDR={h.p_val_adj:.2e}"
            for h in hits.itertuples(index=False)
        )
        evidence = "Strong: gene is a marker of an LDSC-enriched cell type"
        if ldsc_fdr is None or ldsc_fdr >= 0.05:
            evidence = "Moderate: gene is a scRNA marker, LDSC support not significant"

        records.append(
            {
                "Gene": row.Gene,
                "Endometriosis_scRNA_marker_cell_type": best["cell_type"],
                "Endometriosis_scRNA_marker_log2FC": float(best["avg_log2FC"]),
                "Endometriosis_scRNA_marker_FDR": float(best["p_val_adj"]),
                "Endometriosis_all_marker_hits": all_hits,
                "Endometriosis_celltype_LDSC_enrichment": ldsc_enrichment,
                "Endometriosis_celltype_LDSC_FDR": ldsc_fdr,
                "Evidence_level": evidence,
            }
        )

    summary = pd.DataFrame(records)
    keep_cols = [
        "Gene",
        "MAGMA_log10P_fibroids",
        "MAGMA_log10P_endometriosis",
        "MAGMA_P_fibroids",
        "MAGMA_P_endometriosis",
    ]
    return shared[keep_cols].merge(summary, on="Gene", how="left")


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    shared = load_shared_genes()
    endo = endometriosis_marker_summary(shared)
    endo.to_csv(OUT_ENDO, index=False)
    endo.to_csv(ASSET_ENDO, index=False)
    print(f"Wrote {OUT_ENDO}")
    print(f"Wrote {ASSET_ENDO}")


if __name__ == "__main__":
    main()
