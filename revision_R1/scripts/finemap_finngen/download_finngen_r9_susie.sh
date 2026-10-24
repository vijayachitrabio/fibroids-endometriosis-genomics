#!/usr/bin/env bash
# Public FinnGen R9 SuSiE fine-mapping summaries (in-sample LD), used as the primary fine-mapping results.
# Not redistributed here; download from the FinnGen public bucket:
set -euo pipefail
B=https://storage.googleapis.com/finngen-public-data-r9/finemapping/summary
for p in CD2_BENIGN_LEIOMYOMA_UTERI N14_ENDOMETRIOSIS; do
  for f in SUSIE.cred.summary.tsv SUSIE.snp.filter.tsv; do
    curl -fsSLO "$B/finngen_R9_$p.$f"
  done
done
# Comparison with the 1000 Genomes-LD SuSiE: see compare_finngen_insample.py
