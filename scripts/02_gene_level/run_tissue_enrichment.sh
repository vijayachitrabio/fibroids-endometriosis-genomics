#!/usr/bin/env bash
# =============================================================================
# run_tissue_enrichment.sh
# MAGMA gene-property analysis — GTEx v8 tissue enrichment
# Uterine Fibroids and Endometriosis
#
# PREREQUISITES (same MAGMA binary used for gene-level analysis):
#   magma              — already in this folder from previous run
#   NCBI37.3.gene.loc  — already in this folder
#
# WHAT THIS DOES:
#   1. Downloads the GTEx v8 tissue expression covariate file (~150 MB)
#      (pre-computed by the MAGMA team; 54 tissues × 19,131 genes)
#   2. Runs MAGMA --gene-property on fibroids + endometriosis .genes.raw files
#   3. Output: results/fibroids_tissue.gsa.out + endometriosis_tissue.gsa.out
#
# RUNTIME:  ~5-10 min per trait
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAGMA_BIN="$SCRIPT_DIR/magma"
RESULTS="$SCRIPT_DIR/results"
GTEX_FILE="$SCRIPT_DIR/gtex_v8_ts_general_avg_log2TPM.tsv"
GTEX_GZ="$GTEX_FILE.gz"

FIB_RAW="$RESULTS/fibroids_genes.genes.raw"
ENDO_RAW="$RESULTS/endometriosis_genes.genes.raw"

echo "============================================================"
echo "  MAGMA Tissue Enrichment — GTEx v8"
echo "  Uterine Fibroids and Endometriosis"
echo "  $(date)"
echo "============================================================"

# ── Pre-flight checks ────────────────────────────────────────
check_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: Missing required file: $1"; exit 1
    fi
    echo "  [OK] $(basename $1)"
}

echo ""
echo "Pre-flight checks..."
check_file "$MAGMA_BIN"
check_file "$FIB_RAW"
check_file "$ENDO_RAW"
chmod +x "$MAGMA_BIN"

# ── Step 1: Download GTEx v8 gene-property covariate file ────
# Source: MAGMA auxiliary files (CTG lab, Vrije Universiteit Amsterdam)
# File: gtex_v8_ts_general_avg_log2TPM.tsv.gz
# Contains: 54 GTEx tissues, genes identified by Entrez ID (GRCh37/hg19)
# log2(TPM+1) averaged across samples per tissue, then z-scored across tissues
echo ""
echo "------------------------------------------------------------"
echo "STEP 1: Download GTEx v8 tissue covariate file"
echo "------------------------------------------------------------"

if [[ -f "$GTEX_FILE" ]]; then
    echo "  [SKIP] GTEx covariate file already present: $GTEX_FILE"
    echo "  Lines: $(wc -l < "$GTEX_FILE") (first line = header)"
else
    echo "  Building GTEx v8 covariate file from official GTEx source..."
    echo "  (Downloads ~50 MB from Google Cloud Storage, takes ~3-5 min)"
    python3 "$SCRIPT_DIR/prepare_gtex_covariate.py"
    if [[ ! -f "$GTEX_FILE" ]]; then
        echo "ERROR: prepare_gtex_covariate.py did not produce $GTEX_FILE"
        exit 1
    fi
fi

# Verify it looks correct: first column should be integers (Entrez IDs)
echo "  File preview (first 2 columns, first 3 rows):"
awk '{print $1, $2}' "$GTEX_FILE" | head -3 | sed 's/^/    /'

# Sanity: first data column should be a number
FIRST_DATA=$(awk 'NR==2{print $1}' "$GTEX_FILE")
if ! [[ "$FIRST_DATA" =~ ^[0-9]+$ ]]; then
    echo "ERROR: GTEx covariate file does not look correct."
    echo "  Expected integer (Entrez ID) in column 1, row 2; found: $FIRST_DATA"
    exit 1
fi
echo "  Format check passed (Entrez ID in col 1: $FIRST_DATA)"

# ── Step 2: MAGMA gene-property analysis ─────────────────────
# --gene-results : .genes.raw from the gene-level analysis step
# --gene-covar   : GTEx tissue expression file (continuous covariate)
# Model: linear regression of gene Z-stat on tissue expression
#        controlling for gene size (NSNPS, NPARAM) and MAC via covariates
#        already embedded in the .genes.raw covariance structure
# Correction: within-analysis Bonferroni (nrow of .gsa.out)
echo ""
echo "------------------------------------------------------------"
echo "STEP 2: MAGMA gene-property (tissue enrichment)"
echo "------------------------------------------------------------"

for TRAIT in fibroids endometriosis; do
    RAW="$RESULTS/${TRAIT}_genes.genes.raw"
    OUT="$RESULTS/${TRAIT}_tissue"
    echo ""
    echo "  Running: $TRAIT"
    "$MAGMA_BIN" \
        --gene-results "$RAW" \
        --gene-covar "$GTEX_FILE" \
        --model direction=pos \
        --out "$OUT"
    echo "  Done: ${OUT}.gsa.out"
done

# ── Step 3: Quick summary ─────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "STEP 3: Summary of significant tissues"
echo "------------------------------------------------------------"
python3 - "$RESULTS" <<'PYEOF'
import sys, os, pandas as pd

results = sys.argv[1]
thresh = 0.05

for trait in ["fibroids", "endometriosis"]:
    f = f"{results}/{trait}_tissue.gsa.out"
    if not os.path.exists(f):
        print(f"  {trait}: output file not found"); continue
    df = pd.read_csv(f, sep=r'\s+', comment='#')
    n_tests = len(df)
    bonf = thresh / n_tests
    sig = df[df['P'] < bonf].sort_values('P')
    print(f"\n  {trait.upper()} — {n_tests} tissues tested (Bonferroni P<{bonf:.4f})")
    if len(sig) > 0:
        for _, row in sig.iterrows():
            print(f"    {row['VARIABLE']:<40}  P={row['P']:.3e}  beta={row['BETA']:.4f}")
    else:
        print("    No tissues significant after Bonferroni correction")
PYEOF

echo ""
echo "============================================================"
echo "  COMPLETE: $(date)"
echo "  Next step: python3 plot_tissue_enrichment.py"
echo "============================================================"
