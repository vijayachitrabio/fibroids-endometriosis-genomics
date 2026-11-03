"""R2 (Reviewer 1, comment 1): LAVA sensitivity tables from the archived per-specification runs.

Writes to revision_reruns/lava_R2/:
  LAVA_all148_rows_by_specification_long.csv   every bivariate row in every specification, with raw P, estimate, CI,
                                               FDR within the run, FDR in a fixed family (BH over the 88 primary tests),
                                               and the reason when no valid P value was returned
  LAVA_stable_core_three_definitions.csv       the 16 primary blocks classified under (i) FDR within each run,
                                               (ii) fixed family, (iii) the 14 blocks estimable in all five runs
  LAVA_common_estimable_set.csv                the blocks with a valid P value in all five specifications
  LAVA_failed_locus_processing.csv             the two blocks that failed in every run
"""
from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests

RR = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns")
L, OUT = RR / "lava", RR / "lava_R2"
OUT.mkdir(exist_ok=True)
RUNS = {  # label -> (results file, qc file)
    "0": ("res_grid_corr_0.0.csv", "qc_grid_corr_0.0.csv"),
    "0.1724": ("LAVA_corrected_N_ov01724_all_bivariate_rows.csv", None),
    "0.259": ("res_grid_corr_0.2590.csv", "qc_grid_corr_0.2590.csv"),
    "0.35": ("res_grid_corr_0.35.csv", "qc_grid_corr_0.35.csv"),
    "0.369": ("res_grid_corr_0.3685.csv", "qc_grid_corr_0.3685.csv"),
}
OVERLAP_AWARE = ["0.1724", "0.259", "0.35", "0.369"]

d = {k: pd.read_csv(L / v[0]) for k, v in RUNS.items()}
valid = {k: set(x.LOC[x.p.notna()]) for k, x in d.items()}
primary = sorted(valid["0.1724"])
assert len(primary) == 88

def bh(p, m=None):
    p = pd.Series(p).reset_index(drop=True)
    m = m or len(p)
    o = p.sort_values().index
    q = pd.Series(index=p.index, dtype=float)
    run = 1.0
    for rank, i in zip(range(len(p), 0, -1), o[::-1]):
        run = min(run, p[i] * m / rank)
        q[i] = run
    return q.values

rows = []
for k, x in d.items():
    x = x.copy()
    v = x[x.p.notna()].copy()
    v["fdr_within_run"] = bh(v.p.values)
    vf = v[v.LOC.isin(primary)].copy()
    vf["fdr_fixed88"] = bh(vf.p.values, m=88)
    # most conservative family: all 148 rows enter, rows without a valid P are assigned P = 1
    a = x[["LOC"]].copy(); a["fdr_fixed148"] = bh(x.p.fillna(1.0).values, m=148)
    x = x.merge(v[["LOC", "fdr_within_run"]], on="LOC", how="left").merge(vf[["LOC", "fdr_fixed88"]], on="LOC", how="left").merge(a, on="LOC", how="left")
    x.loc[x.p.isna(), "fdr_fixed148"] = float("nan")
    x["specification"] = k
    x["valid_P"] = x.p.notna()
    x["reason_no_valid_P"] = ["" if ok else "bivariate estimate returned as NA (outside |rho| <= 1.25 bound or not estimable; LAVA does not report which)"
                              for ok in x.valid_P]
    rows.append(x)
long = pd.concat(rows, ignore_index=True)
long = long[["specification", "LOC", "CHR", "START", "STOP", "n_snps", "K", "univ_p_endo", "univ_p_uf", "h2_endo", "h2_uf",
             "rho", "rho.lower", "rho.upper", "p", "valid_P", "fdr_within_run", "fdr_fixed88", "fdr_fixed148", "reason_no_valid_P"]]
long.to_csv(OUT / "LAVA_all148_rows_by_specification_long.csv", index=False)

fail = pd.read_csv(L / "qc_grid_corr_0.0.csv")
fail = fail[fail.status != "ok"]
fail.to_csv(OUT / "LAVA_failed_locus_processing.csv", index=False)

common = sorted(set.intersection(*valid.values()))
pd.DataFrame({"LOC": common}).to_csv(OUT / "LAVA_common_estimable_set.csv", index=False)

# classification of the 16 primary blocks
w = {k: long[long.specification == k].set_index("LOC") for k in RUNS}
sig16 = sorted(l for l in primary if w["0.1724"].loc[l, "fdr_within_run"] < 0.05)
assert len(sig16) == 16
cm = {}
for k in RUNS:
    sub = w[k].loc[common]
    cm[k] = pd.Series(bh(sub.p.values), index=common)
out = []
for l in sig16:
    r = {"Block": l, "CHR": int(w["0.1724"].loc[l, "CHR"]), "START": int(w["0.1724"].loc[l, "START"]), "STOP": int(w["0.1724"].loc[l, "STOP"])}
    for k in RUNS:
        r[f"estimable [{k}]"] = l in valid[k]
        r[f"FDR within run [{k}]"] = w[k].loc[l, "fdr_within_run"] if l in valid[k] else float("nan")
        r[f"FDR fixed family m=88 [{k}]"] = w[k].loc[l, "fdr_fixed88"] if l in valid[k] else float("nan")
        r[f"FDR common set [{k}]"] = cm[k].get(l, float("nan"))
        r[f"FDR fixed family m=148 [{k}]"] = w[k].loc[l, "fdr_fixed148"] if l in valid[k] else float("nan")
    def stable(col):
        ok = all(r[f"{col} [{k}]"] < 0.05 for k in OVERLAP_AWARE)
        z = r[f"{col} [0]"]
        return ok and (z != z or z < 0.05)   # nan = not estimable without overlap
    r["stable: within-run FDR"] = stable("FDR within run")
    r["stable: fixed family m=148"] = stable("FDR fixed family m=148")
    r["stable: fixed family"] = stable("FDR fixed family m=88")
    r["tested without overlap"] = l in valid["0"]
    r["in common estimable set"] = l in common
    r["stable: strict (tested at 0 and significant everywhere, fixed family)"] = r["stable: fixed family"] and r["tested without overlap"]
    out.append(r)
o = pd.DataFrame(out)
o.to_csv(OUT / "LAVA_stable_core_three_definitions.csv", index=False)

pd.set_option("display.width", 250)
print("valid tests:", {k: len(v) for k, v in valid.items()}, "common:", len(common))
print(o[["Block", "CHR", "stable: within-run FDR", "stable: fixed family", "tested without overlap", "in common estimable set",
         "stable: strict (tested at 0 and significant everywhere, fixed family)"]].to_string())
for col in ("stable: within-run FDR", "stable: fixed family", "stable: fixed family m=148", "stable: strict (tested at 0 and significant everywhere, fixed family)"):
    print(col, int(o[col].sum()))
print("common set: primary-16 members:", int(o["in common estimable set"].sum()))
for k in RUNS:
    print(k, "sig within common set:", int((cm[k] < 0.05).sum()), sorted(cm[k][cm[k] < 0.05].index))
