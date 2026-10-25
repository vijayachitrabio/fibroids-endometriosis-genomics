#!/usr/bin/env bash
# S-LDSC with an all-genes control annotation (Finucane et al. 2018): every GSE218044 gene with a
# GRCh37 coordinate (NCBI37.3; Ensembl for marker genes absent from it), +/-100 kb, fitted jointly with
# the baseline model and each cell-type annotation. Marker sets: original and confident-cell.
set -uo pipefail
R=/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns/sldsc_reannot
ECLAI=/Users/vijayachitramodhukur/Desktop/ECLAI
LDSC="python3 $R/ldsc_tool/ldsc.py"
BED=$ECLAI/pregnancy_chronic_continuity_audit/stage3/reference/1000G.EUR.QC.
REFS=$ECLAI/uterine_fibroids/ldsc_fibroids_rerun/refs
W=$ECLAI/scEndoExplorer_v1/ldsc_refs/ldsc_refs/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=$R/frq/1000G.EUR.QC.
CTRL=$R/annots_control
gwas() { case $1 in fibroids) echo $ECLAI/uterine_fibroids/outputs_v2/Uterine_Fibroids_munged.sumstats.gz ;;
                     endometriosis) echo $ECLAI/uterine_fibroids/outputs_v2/Endometriosis_munged.sumstats.gz ;; esac; }
CTS="B_cell Endothelial Epithelial Macrophage Smooth_muscle Stromal T_cell__CD8_"
for C in $(seq 1 22); do
  [ -f $CTRL/AllGenes.$C.l2.ldscore.gz ] && continue
  python3 $R/scripts/make_annot.py $CTRL/AllGenes.genes $CTRL/AllGenes_coord.txt ${BED}$C.bim $CTRL/AllGenes.$C.annot.gz > /dev/null
  $LDSC --l2 --bfile ${BED}$C --ld-wind-cm 1 --thin-annot --annot $CTRL/AllGenes.$C.annot.gz \
    --print-snps $R/hm3_snps.txt --out $CTRL/AllGenes.$C > $CTRL/AllGenes.$C.l2.log 2>&1
done
echo "control LD scores done"
for SET in original confident; do
  ANN=$([ $SET = original ] && echo $ECLAI/uterine_fibroids/ldsc_fibroids_rerun/annots || echo $R/annots_confident)
  for T in fibroids endometriosis; do
    for M in v12 ld22; do
      B=$([ $M = v12 ] && echo $REFS/baseline_v1.2/baseline. || echo $REFS/baselineLD_v2.2/baselineLD.)
      O=$R/results_control_${SET}_${T}_$M; mkdir -p $O
      for CT in $CTS; do
        $LDSC --h2 $(gwas $T) --ref-ld-chr ${B}@,$CTRL/AllGenes.@,$ANN/$CT.@ --w-ld-chr ${W}@ --frqfile-chr ${FRQ}@ \
          --overlap-annot --print-coefficients --out $O/$CT > /dev/null 2>&1
      done
      echo "control $SET $T $M done"
    done
  done
done
echo "CONTROL DONE"
