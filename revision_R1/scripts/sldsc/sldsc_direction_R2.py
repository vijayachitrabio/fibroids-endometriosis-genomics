"""R2 (Reviewer 1, comment 2): dependence of the S-LDSC stromal result on marker assignment and test direction.

One-sided coefficient P comes from the S-LDSC output (coefficient Z, upper tail); two-sided P = 2 * min(p, 1 - p).
FDR: Benjamini-Hochberg over the seven annotations within trait, model and marker set.
Writes revision_reruns/sldsc_R2/SLDSC_direction_and_marker_dependence.csv and prints the stromal rows.
"""
from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests

RR = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns")
OUT = RR / "sldsc_R2"
OUT.mkdir(exist_ok=True)
d = pd.read_csv(RR / "sldsc_reannot/SLDSC_original_vs_confident_markers.csv")
d["coef_p_two"] = 2 * d.coef_p.clip(upper=1 - d.coef_p)
g = ["set", "trait", "model"]
d["coef_fdr_one_sided"] = d.groupby(g).coef_p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
d["coef_fdr_two_sided"] = d.groupby(g).coef_p_two.transform(lambda p: multipletests(p, method="fdr_bh")[1])
d["set"] = d.set.map({"original": "original markers", "confident": "re-derived markers (concordant cells)"}).fillna(d.set)
d.to_csv(OUT / "SLDSC_direction_and_marker_dependence.csv", index=False)
pd.set_option("display.width", 250)
print(d[d.ct.str.contains("Strom")][["set", "trait", "model", "prop", "enr", "enr_fdr", "z", "coef_p", "coef_fdr", "coef_fdr_one_sided", "coef_fdr_two_sided"]].round(4).to_string())
print(sorted(d.ct.unique()), d.set.unique())
