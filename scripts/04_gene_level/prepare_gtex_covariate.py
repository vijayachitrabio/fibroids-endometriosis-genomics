#!/usr/bin/env python3
"""
prepare_gtex_covariate.py
=========================
Build the GTEx v8 tissue covariate file required by MAGMA gene-property analysis.

SOURCE:  Official GTEx v8 median-TPM bulk download (Google Cloud Storage, public)
         GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz  (~50 MB)

OUTPUT:  gtex_v8_ts_general_avg_log2TPM.tsv  (same name/format MAGMA expects)
         - Tab-separated, one header row
         - Column 1: Entrez gene ID (matches MAGMA .genes.raw)
         - Columns 2-55: log2(median_TPM + 1) per GTEx tissue
         - 54 tissues, ~17,000–18,000 genes (those mappable to Entrez IDs)

GENE ID MAPPING:
         Entrez IDs come from NCBI37.3.gene.loc (already in this folder).
         GTEx uses ENSEMBL IDs + gene symbols; we join on gene symbol.
         Duplicate symbols → keep first occurrence (rare, <0.5 % of genes).

USAGE:   python3 prepare_gtex_covariate.py
         Run from the magma_inputs/ directory.

RUNTIME: ~3-5 min (dominated by 50 MB download + decompression)
"""

import os
import sys
import gzip
import shutil
import urllib.request

import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GENE_LOC   = os.path.join(SCRIPT_DIR, "NCBI37.3.gene.loc")
GCT_GZ     = os.path.join(SCRIPT_DIR, "GTEx_v8_median_tpm.gct.gz")
GCT_FILE   = os.path.join(SCRIPT_DIR, "GTEx_v8_median_tpm.gct")
OUT_FILE   = os.path.join(SCRIPT_DIR, "gtex_v8_ts_general_avg_log2TPM.tsv")

GTEX_URL = (
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)


# ── Step 1: Download GTEx median TPM ─────────────────────────────────────────
def download_with_progress(url, dest):
    """Download url to dest, showing progress."""
    def reporthook(count, block_size, total_size):
        if total_size > 0:
            pct = min(int(count * block_size * 100 / total_size), 100)
            sys.stdout.write(f"\r  Progress: {pct}%  ")
            sys.stdout.flush()
    print(f"  URL: {url}")
    urllib.request.urlretrieve(url, dest, reporthook)
    print()  # newline after progress


print("=" * 60)
print("  GTEx v8 covariate file builder for MAGMA")
print("=" * 60)

# Check NCBI37.3.gene.loc exists
if not os.path.exists(GENE_LOC):
    sys.exit(f"ERROR: {GENE_LOC} not found.\n"
             f"Please ensure NCBI37.3.gene.loc is in the magma_inputs/ folder.")

print(f"\nStep 1: Download GTEx v8 median TPM")
if os.path.exists(GCT_FILE):
    print(f"  [SKIP] Uncompressed GCT already present: {GCT_FILE}")
elif os.path.exists(GCT_GZ) and os.path.getsize(GCT_GZ) > 1_000_000:
    print(f"  [SKIP] GCT gz already downloaded: {GCT_GZ}")
else:
    print(f"  Downloading (~50 MB) ...")
    download_with_progress(GTEX_URL, GCT_GZ)
    sz = os.path.getsize(GCT_GZ)
    print(f"  Downloaded: {sz / 1e6:.1f} MB")
    if sz < 1_000_000:
        sys.exit("ERROR: Downloaded file too small — download may have failed.\n"
                 f"  File: {GCT_GZ}\n"
                 f"  Check your internet connection and retry.")

if not os.path.exists(GCT_FILE):
    print(f"  Decompressing ...")
    with gzip.open(GCT_GZ, 'rb') as f_in, open(GCT_FILE, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"  Done: {GCT_FILE}")


# ── Step 2: Parse GTEx GCT file ───────────────────────────────────────────────
# GCT format:
#   Line 1: #1.2
#   Line 2: <n_genes>  <n_samples>
#   Line 3: Name  Description  <tissue1>  <tissue2>  ...   (header)
#   Lines 4+: ENSG00000xxxxx.N  SYMBOL  <tpm1>  <tpm2>  ...
print(f"\nStep 2: Parse GTEx GCT file")
with open(GCT_FILE, 'r') as f:
    _line1 = f.readline()          # #1.2
    _line2 = f.readline()          # n_genes  n_tissues
    header = f.readline().rstrip('\n').split('\t')

tissue_cols = header[2:]           # skip Name, Description
print(f"  Tissues found: {len(tissue_cols)}")
print(f"  Sample tissue names: {tissue_cols[:4]} ...")

# Read data (skip 3-line GCT header)
print(f"  Reading expression matrix (this takes ~30 sec) ...")
gtex = pd.read_csv(GCT_FILE, sep='\t', skiprows=2,
                   usecols=['Name', 'Description'] + tissue_cols,
                   low_memory=False)
gtex.rename(columns={'Name': 'ENSEMBL', 'Description': 'SYMBOL'}, inplace=True)
# Strip version from ENSEMBL ID (ENSG00000xxx.14 → ENSG00000xxx)
gtex['ENSEMBL'] = gtex['ENSEMBL'].str.split('.').str[0]
print(f"  Loaded: {len(gtex):,} genes × {len(tissue_cols)} tissues")


# ── Step 3: Load Entrez↔Symbol mapping from NCBI37.3.gene.loc ────────────────
# Format: ENTREZ_ID  CHR  START  STOP  STRAND  SYMBOL
print(f"\nStep 3: Load Entrez ID mapping from {os.path.basename(GENE_LOC)}")
gene_loc = pd.read_csv(
    GENE_LOC, sep='\t', header=None,
    names=['ENTREZ', 'CHR', 'START', 'STOP', 'STRAND', 'SYMBOL'],
    usecols=['ENTREZ', 'SYMBOL']
)
gene_loc.dropna(subset=['SYMBOL'], inplace=True)
gene_loc['SYMBOL'] = gene_loc['SYMBOL'].str.strip()
gene_loc.drop_duplicates(subset='SYMBOL', keep='first', inplace=True)
print(f"  Unique gene symbols in NCBI37.3.gene.loc: {len(gene_loc):,}")


# ── Step 4: Merge on gene symbol ─────────────────────────────────────────────
print(f"\nStep 4: Merge GTEx expression with Entrez IDs on gene symbol")
gtex['SYMBOL'] = gtex['SYMBOL'].str.strip()
merged = gtex.merge(gene_loc, on='SYMBOL', how='inner')
print(f"  Matched genes: {len(merged):,} "
      f"(of {len(gtex):,} GTEx × {len(gene_loc):,} NCBI37)")
print(f"  Unmatched (not in NCBI37.3): {len(gtex) - len(merged):,} "
      f"(mostly novel genes, pseudogenes, non-coding — expected)")

if len(merged) < 10_000:
    print("WARNING: Unexpectedly few matches — check gene symbol format.")


# ── Step 5: log2(TPM + 1) transform ──────────────────────────────────────────
print(f"\nStep 5: log2(TPM + 1) transform")
expr_cols = tissue_cols
merged[expr_cols] = np.log2(merged[expr_cols].values.astype(float) + 1.0)

# Verify no NaN or infinite values
n_nan = merged[expr_cols].isnull().sum().sum()
n_inf = np.isinf(merged[expr_cols].values).sum()
print(f"  NaN values: {n_nan}  |  Inf values: {n_inf}")
if n_nan > 0 or n_inf > 0:
    merged[expr_cols] = merged[expr_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    print("  Fixed: replaced NaN/Inf with 0.0")


# ── Step 6: Write MAGMA-compatible covariate file ─────────────────────────────
# Format: GENE<TAB>Tissue1<TAB>Tissue2...  (Entrez IDs, no ENSG/symbol)
print(f"\nStep 6: Write output file")
out_df = merged[['ENTREZ'] + expr_cols].copy()
out_df.rename(columns={'ENTREZ': 'GENE'}, inplace=True)
out_df = out_df.drop_duplicates(subset='GENE', keep='first')
out_df = out_df.sort_values('GENE')
out_df.to_csv(OUT_FILE, sep='\t', index=False, float_format='%.6f')

sz = os.path.getsize(OUT_FILE)
print(f"  Written: {OUT_FILE}")
print(f"  Rows: {len(out_df):,} genes  |  Columns: {out_df.shape[1]}  "
      f"(1 ID + {out_df.shape[1]-1} tissues)  |  Size: {sz/1e6:.1f} MB")

# Sanity check: preview
print(f"\n  First 3 rows (GENE + first 3 tissues):")
print(out_df.iloc[:3, :4].to_string(index=False))

# Clean up large intermediate GCT file to save disk space
print(f"\nStep 7: Cleanup intermediate files")
if os.path.exists(GCT_FILE):
    os.remove(GCT_FILE)
    print(f"  Removed: {GCT_FILE}")
if os.path.exists(GCT_GZ):
    os.remove(GCT_GZ)
    print(f"  Removed: {GCT_GZ}")

print(f"\n{'=' * 60}")
print(f"  Done. GTEx covariate file ready for MAGMA:")
print(f"  {OUT_FILE}")
print(f"  Next step: bash run_tissue_enrichment.sh")
print(f"{'=' * 60}")
