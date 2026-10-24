#!/usr/bin/env bash
# LDSC genetic correlation under three sample-size inputs (Reviewer 1, comment 1).
# The munged FinnGen R9 files carry a constant N; only the N column is replaced,
# Z scores are unchanged. rg is invariant to a constant rescaling of N.
#   submitted : constant N used in the original analysis (218,728 / 239,335)
#   total     : R9 endpoint totals (210,870 / 122,652)
#   neff      : 4 / (1/cases + 1/controls) (107,629 / 52,928)
# Inputs: Uterine_Fibroids_munged.sumstats.gz, Endometriosis_munged.sumstats.gz (LDSC-munged
# FinnGen R9 CD2_BENIGN_LEIOMYOMA_UTERI and N14_ENDOMETRIOSIS), eur_w_ld_chr/ LD scores.
# Output logs: results/ldsc/rg_{submitted,total,neff}.log
set -euo pipefail
MUNGED=${MUNGED:-outputs_v2}; LDSC=${LDSC:-ldsc.py}; W=${W:-eur_w_ld_chr/}
python3 - <<PY
import pandas as pd
for k, f, c, n in (("UF", "$MUNGED/Uterine_Fibroids_munged.sumstats.gz", 31661, 179209),
                   ("ENDO", "$MUNGED/Endometriosis_munged.sumstats.gz", 15088, 107564)):
    d = pd.read_csv(f, sep="\t")
    d.to_csv(f"{k}_submitted.sumstats.gz", sep="\t", index=False, compression="gzip")
    for tag, N in (("total", c + n), ("neff", 4 / (1 / c + 1 / n))):
        e = d.copy(); e["N"] = N
        e.to_csv(f"{k}_{tag}.sumstats.gz", sep="\t", index=False, compression="gzip")
PY
for s in submitted total neff; do
  python3 "$LDSC" --rg ENDO_$s.sumstats.gz,UF_$s.sumstats.gz --ref-ld-chr "$W" --w-ld-chr "$W" --out rg_$s
done
