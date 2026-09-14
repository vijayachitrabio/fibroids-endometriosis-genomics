#!/usr/bin/env bash
# =============================================================================
# run_magma.sh — Full MAGMA gene-level analysis for fibroids & endometriosis
# =============================================================================
#
# PREREQUISITES (copy these files to this directory before running):
#   magma              — MAGMA binary  (download from ctg.cncr.nl/software/magma)
#   NCBI37.3.gene.loc  — Gene location file (bundled with MAGMA download)
#
# FILES ALREADY PREPARED IN THIS DIRECTORY:
#   fibroids.snploc        — SNP  CHR  BP  (12,018,721 SNPs)
#   fibroids.pval          — SNP  P    N   (12,018,721 SNPs)
#   endometriosis.snploc   — SNP  CHR  BP  (12,015,108 SNPs)
#   endometriosis.pval     — SNP  P    N   (12,015,108 SNPs)
#
# LD REFERENCE (already on disk):
#   ../mr_analysis_2026_04_16/reference_data/g1000_eur  (.bed/.bim/.fam)
#   22,665,064 SNPs, 503 EUR samples (1000G Phase 3)
#
# OUTPUT:
#   results/fibroids_genes.genes.out
#   results/endometriosis_genes.genes.out
#   results/fibroids_vs_endo_comparison.tsv   (merged results)
#
# RUNTIME ESTIMATE:  ~45-90 min per trait on a modern laptop
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LD_REF="$SCRIPT_DIR/../mr_analysis_2026_04_16/reference_data/g1000_eur"
GENE_LOC="$SCRIPT_DIR/NCBI37.3.gene.loc"
MAGMA_BIN="$SCRIPT_DIR/magma"
RESULTS="$SCRIPT_DIR/results"
ANNOT="$SCRIPT_DIR/annot"

latest_match() {
    python3 - "$SCRIPT_DIR" "$1" "$2" <<'PYEOF'
import glob
import os
import sys

base = sys.argv[1]
label = sys.argv[2]
suffix = sys.argv[3]
patterns = {
    ("fibroids", "snploc"): ["fibroids.snploc", "*fibroid*.snploc"],
    ("fibroids", "pval"): ["fibroids.pval", "*fibroid*.pval"],
    ("endometriosis", "snploc"): ["endometriosis.snploc", "*endometr*.snploc"],
    ("endometriosis", "pval"): ["endometriosis.pval", "*endometr*.pval"],
}

candidates = []
for pattern in patterns[(label, suffix)]:
    candidates.extend(glob.glob(os.path.join(base, pattern)))

candidates = [p for p in candidates if os.path.isfile(p)]
candidates = sorted(set(candidates), key=lambda p: (os.path.getmtime(p), p), reverse=True)
print(candidates[0] if candidates else "")
PYEOF
}

FIB_SNPloc="$(latest_match fibroids snploc)"
FIB_PVAL="$(latest_match fibroids pval)"
ENDO_SNPloc="$(latest_match endometriosis snploc)"
ENDO_PVAL="$(latest_match endometriosis pval)"

# ─────────────────────────────────────────────────────────────
# Pre-flight checks
# ─────────────────────────────────────────────────────────────
echo "============================================================"
echo "  MAGMA gene-level analysis — fibroids & endometriosis"
echo "  $(date)"
echo "============================================================"

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: Required file missing: $1"
        echo "  Please read the header of this script for setup instructions."
        exit 1
    fi
    echo "  [OK] $1"
}

echo ""
echo "Pre-flight checks..."
check_file "$MAGMA_BIN"
check_file "$GENE_LOC"
check_file "$LD_REF.bed"
check_file "$LD_REF.bim"
check_file "$LD_REF.fam"
check_file "$FIB_SNPloc"
check_file "$FIB_PVAL"
check_file "$ENDO_SNPloc"
check_file "$ENDO_PVAL"

chmod +x "$MAGMA_BIN"
mkdir -p "$RESULTS" "$ANNOT"
echo "All checks passed."
echo ""
echo "Using latest MAGMA inputs:"
echo "  fibroids.snploc      -> $(basename "$FIB_SNPloc")"
echo "  fibroids.pval        -> $(basename "$FIB_PVAL")"
echo "  endometriosis.snploc -> $(basename "$ENDO_SNPloc")"
echo "  endometriosis.pval   -> $(basename "$ENDO_PVAL")"
echo ""

# ─────────────────────────────────────────────────────────────
# STEP 1: SNP annotation  (map SNPs → genes using window ±0 bp)
# ─────────────────────────────────────────────────────────────
echo "------------------------------------------------------------"
echo "STEP 1: SNP annotation (gene window 0 kb)"
echo "------------------------------------------------------------"

for TRAIT in fibroids endometriosis; do
    if [[ "$TRAIT" == "fibroids" ]]; then
        SNPLOC="$FIB_SNPloc"
    else
        SNPLOC="$ENDO_SNPloc"
    fi
    echo ""
    echo "  Annotating $TRAIT ..."
    "$MAGMA_BIN" \
        --annotate \
        --snp-loc "$SNPLOC" \
        --gene-loc "$GENE_LOC" \
        --out "$ANNOT/${TRAIT}"
    echo "  Done: $ANNOT/${TRAIT}.genes.annot"
done

# ─────────────────────────────────────────────────────────────
# STEP 2: Gene-level analysis  (SNP p-values → gene p-values)
# Using mean model (default) + 1000G EUR LD reference
# ─────────────────────────────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "STEP 2: Gene-level analysis"
echo "------------------------------------------------------------"

for TRAIT in fibroids endometriosis; do
    if [[ "$TRAIT" == "fibroids" ]]; then
        PVAL="$FIB_PVAL"
    else
        PVAL="$ENDO_PVAL"
    fi
    N=$(awk 'NR==2{print $3}' "$PVAL")
    echo ""
    echo "  Running $TRAIT  (N=$N) ..."
    "$MAGMA_BIN" \
        --bfile "$LD_REF" \
        --pval "$PVAL" use=SNP,P ncol=N \
        --gene-annot "$ANNOT/${TRAIT}.genes.annot" \
        --out "$RESULTS/${TRAIT}_genes"
    echo "  Done: $RESULTS/${TRAIT}_genes.genes.out"
done

# ─────────────────────────────────────────────────────────────
# STEP 3: Quick summary report
# ─────────────────────────────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "STEP 3: Summary"
echo "------------------------------------------------------------"
python3 - "$RESULTS" <<'PYEOF'
import os
import sys

import pandas as pd

results_dir = sys.argv[1]
base_dir = os.path.dirname(results_dir)
gene_loc = os.path.join(base_dir, "NCBI37.3.gene.loc")

gene_map = pd.read_csv(
    gene_loc,
    sep="\t",
    header=None,
    names=["GENE", "CHR_LOC", "START_LOC", "STOP_LOC", "STRAND", "SYMBOL"],
)
gene_map = gene_map[["GENE", "SYMBOL"]].drop_duplicates()

TRAITS = ["fibroids", "endometriosis"]
dfs = {}
for trait in TRAITS:
    f = f"{results_dir}/{trait}_genes.genes.out"
    if not os.path.exists(f):
        print(f"  {trait}: results file not found ({f})")
        continue
    df = pd.read_csv(f, sep=r'\s+', comment="#")
    df = df.merge(gene_map, on="GENE", how="left")
    df["TRAIT"] = trait
    dfs[trait] = df
    sig = (df["P"] < 2.5e-6).sum()   # Bonferroni ~0.05/19,993 genes
    print(f"  {trait}: {len(df):,} genes tested, {sig} significant (P<2.5e-6)")

if len(dfs) == 2:
    merged = dfs["fibroids"][["GENE","SYMBOL","P"]].merge(
        dfs["endometriosis"][["GENE","SYMBOL","P"]],
        on=["GENE","SYMBOL"], suffixes=("_fib","_endo")
    )
    merged["shared_sig"] = (merged["P_fib"] < 2.5e-6) & (merged["P_endo"] < 2.5e-6)
    out = f"{results_dir}/fibroids_vs_endo_comparison.tsv"
    merged.to_csv(out, sep="\t", index=False)
    print(f"\n  Shared significant genes: {merged['shared_sig'].sum()}")
    print(f"  Comparison table: {out}")
PYEOF

echo ""
echo "============================================================"
echo "  MAGMA analysis complete:  $(date)"
echo "  Results in:  $RESULTS/"
echo "============================================================"
