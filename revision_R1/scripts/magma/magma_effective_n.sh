#!/usr/bin/env bash
# MAGMA v1.10 gene analysis (SNP-wise mean) with R9 effective sample sizes (Reviewer 1, comment 1).
# Inputs: fibroids.pval / endometriosis.pval (SNP, P from FinnGen R9), gene annotation files
# from MAGMA --annotate (NCBI37.3, 0 kb window), 1000 Genomes Phase 3 EUR PLINK reference (g1000_eur).
# Output: results/magma/{fibroids,endo}_neff.genes.out
set -euo pipefail
MAGMA=${MAGMA:-magma}; REF=${REF:-g1000_eur}
$MAGMA --bfile $REF --pval fibroids.pval use=SNP,P N=107629 --gene-annot annot/fibroids.genes.annot --out fibroids_neff
$MAGMA --bfile $REF --pval endometriosis.pval use=SNP,P N=52928 --gene-annot annot/endometriosis.genes.annot --out endo_neff
