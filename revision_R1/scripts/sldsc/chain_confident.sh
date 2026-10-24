#!/bin/bash
# wait for the marker re-derivation, write gene lists, then run S-LDSC on the confident-cell markers
R=/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns
until grep -q "| done" $R/reannotation/reannotate.log; do
  pgrep -f reannotate_GSE218044 >/dev/null || { echo "reannotation stopped without finishing"; exit 1; }
  sleep 30
done
mkdir -p $R/sldsc_reannot/annots_confident
python3 - <<PY
import pandas as pd
m = pd.read_csv("$R/reannotation/markers_confident_cells.csv")
for ct, g in m.groupby("cell_type"):
    safe = ct.replace(" ", "_").replace("(", "_").replace(")", "_")
    open(f"$R/sldsc_reannot/annots_confident/{safe}.genes", "w").write("\n".join(g.gene.unique()))
    print(safe, g.gene.nunique())
PY
bash $R/sldsc_reannot/scripts/run_sldsc.sh confident
echo "CHAIN DONE"
