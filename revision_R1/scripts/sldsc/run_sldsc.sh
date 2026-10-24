#!/usr/bin/env bash
# S-LDSC for the marker-set sensitivity analysis (Reviewer 1, comment 4).
#   SET=original   existing annotations (validation: must reproduce published results)
#   SET=confident  markers re-derived from cells whose per-cell label matches the cluster label
# Usage: run_sldsc.sh SET
set -uo pipefail
SET=$1
R=/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns/sldsc_reannot
ECLAI=/Users/vijayachitramodhukur/Desktop/ECLAI
LDSC="python3 $R/ldsc_tool/ldsc.py"
BED=$ECLAI/pregnancy_chronic_continuity_audit/stage3/reference/1000G.EUR.QC.
REFS=$ECLAI/uterine_fibroids/ldsc_fibroids_rerun/refs
V12=$REFS/baseline_v1.2/baseline.
LD22=$REFS/baselineLD_v2.2/baselineLD.
W=$ECLAI/scEndoExplorer_v1/ldsc_refs/ldsc_refs/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=$R/frq/1000G.EUR.QC.
gwas() { case $1 in fibroids) echo $ECLAI/uterine_fibroids/outputs_v2/Uterine_Fibroids_munged.sumstats.gz ;;
                     endometriosis) echo $ECLAI/uterine_fibroids/outputs_v2/Endometriosis_munged.sumstats.gz ;; esac; }
CTS="B_cell Endothelial Epithelial Macrophage Smooth_muscle Stromal T_cell__CD8_"

if [ "$SET" = "original" ]; then
  ANN=$ECLAI/uterine_fibroids/ldsc_fibroids_rerun/annots
else
  ANN=$R/annots_$SET
  mkdir -p $ANN
  for CT in $CTS; do
    for C in $(seq 1 22); do
      [ -f $ANN/$CT.$C.l2.ldscore.gz ] && continue
      python3 $R/scripts/make_annot.py $ANN/$CT.genes $ECLAI/scEndoExplorer_v1/ldsc_refs/ldsc_refs/ENSG_coord.txt \
        ${BED}$C.bim $ANN/$CT.$C.annot.gz > /dev/null
      $LDSC --l2 --bfile ${BED}$C --ld-wind-cm 1 --thin-annot --annot $ANN/$CT.$C.annot.gz \
        --print-snps $R/hm3_snps.txt --out $ANN/$CT.$C > $ANN/$CT.$C.l2.log 2>&1
    done
    echo "$CT LD scores done"
  done
fi
for T in fibroids endometriosis; do
  for M in v12 ld22; do
    B=$([ $M = v12 ] && echo $V12 || echo $LD22)
    O=$R/results_${SET}_${T}_$M; mkdir -p $O
    for CT in $CTS; do
      $LDSC --h2 $(gwas $T) --ref-ld-chr ${B}@,$ANN/$CT.@ --w-ld-chr ${W}@ --frqfile-chr ${FRQ}@ \
        --overlap-annot --print-coefficients --out $O/$CT > /dev/null 2>&1
    done
    echo "$SET $T $M done"
  done
done
