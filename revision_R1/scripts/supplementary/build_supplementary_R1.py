"""Build the revised supplementary workbook (R1).

Output: Funciton_int_ge/Supplementary_Tables_S1_to_S6_F_IG_R1.xlsx (new file;
the submitted workbook is not modified).

  S1  Data sources, sample sizes, genome-wide and primary local correlation
  S2  MAGMA gene-level results (relabelled; effective-N rerun; cross-trait)
  S3  Fine-mapping (was S4): FinnGen in-sample SuSiE, 1000G sensitivity
  S4  S-LDSC and single-cell context (was S5)
  S5  LAVA sensitivity across sampling-correlation values (new)
  S6  Genome-build traceability (new)

Former S3 (Reactome terms) is removed.
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

FIG = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge")
RR = FIG / "revision_reruns"
UF = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids")
OLD = FIG / "Supplementary_Tables_S1_to_S5_F_IG.xlsx"
OUT = FIG / "Supplementary_Tables_S1_to_S6_F_IG_R1.xlsx"

UF_CASES, UF_CTRL = 31661, 179209
EN_CASES, EN_CTRL = 15088, 107564
neff = lambda a, b: 4 / (1 / a + 1 / b)

sheets = {}          # name -> list of (panel title, DataFrame)


def panel(sheet, title, df):
    sheets.setdefault(sheet, []).append((title, df.reset_index(drop=True)))


# ======================================================================= S1
prov = pd.DataFrame([
    ["Uterine fibroids", "CD2_BENIGN_LEIOMYOMA_UTERI", "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz",
     UF_CASES, UF_CTRL, UF_CASES + UF_CTRL, round(neff(UF_CASES, UF_CTRL), 1),
     "35,474 / 183,254 (total 218,728)", "LDSC and MAGMA (constant N in munged files)",
     "LAVA: R9 cases/controls; SuSiE: R9 effective N; LDSC and MAGMA: munged-file N (identical results with R9 effective N, Panel B and Table S2)"],
    ["Endometriosis", "N14_ENDOMETRIOSIS", "finngen_R9_N14_ENDOMETRIOSIS.gz",
     EN_CASES, EN_CTRL, EN_CASES + EN_CTRL, round(neff(EN_CASES, EN_CTRL), 1),
     "17,045 / 222,290 (total 239,335)", "LDSC and MAGMA (constant N in munged files)",
     "LAVA: R9 cases/controls; SuSiE: R9 effective N; LDSC and MAGMA: munged-file N (identical results with R9 effective N, Panel B and Table S2)"],
], columns=["Trait", "FinnGen R9 endpoint", "Summary-statistic file", "Cases (R9 endpoint)",
            "Controls (R9 endpoint)", "Total N", "Effective N = 4/(1/cases+1/controls)",
            "Counts carried in the munged input files", "Where used",
            "Sample sizes used in this study"])
panel("S1", "Panel A. Data sources and sample-size provenance. FinnGen R9 files do not report per-variant N; counts are from the R9 endpoint metadata.", prov)

rg = []
for tag, lab in (("submitted", "Munged-file constant N"), ("total", "R9 total N"), ("neff", "R9 effective N")):
    t = open(RR / f"ldsc/rg_{tag}.log").read().splitlines()
    i = [k for k, l in enumerate(t) if l.startswith("Summary of Genetic Correlation")][0]
    v = t[i + 2].split()
    rg.append([lab, "standard LDSC (eur_w_ld_chr)", float(v[2]), float(v[3]), float(v[4]), float(v[5]),
               float(v[6]), float(v[10])])
rg.append(["Munged-file constant N (estimate reported in the main text)", "GenomicSEM multivariable LDSC", 0.511, 0.0691,
           np.nan, 1.37e-13, np.nan, 0.1917])
panel("S1", "Panel B. Genome-wide genetic correlation (LDSC) under alternative sample sizes. rg does not depend on a constant rescaling of N; observed-scale h2 does.",
      pd.DataFrame(rg, columns=["Sample-size input", "Software", "rg", "SE", "Z", "P",
                                "Observed-scale h2 (fibroids)", "Cross-trait intercept"]))

elig = pd.DataFrame([
    ["Blocks analysed", 2495, 2495],
    ["Failed during locus processing", 2, 2],
    ["No bivariate estimate (local h2 not estimable for one or both traits)", 2343, 2345],
    ["Reached the bivariate test", 150, 148],
    ["Estimate outside parameter bound (|rho| > 1.25), P set to missing", 58, 60],
    ["Valid bivariate tests", 92, 88],
    ["FDR-significant blocks (Benjamini-Hochberg over valid tests)", 16, 16],
], columns=["Step", "Munged-file case/control counts", "R9 case/control counts (primary)"])
panel("S1", "Panel C. LAVA eligibility (2,495 LD blocks, 1000 Genomes EUR reference, sampling correlation 0.1724, univ.thresh = 1)", elig)

lv = pd.read_csv(RR / "lava/LAVA_corrected_N_ov01724_all_valid.csv")
lv = lv[["LOC", "CHR", "START", "STOP", "n_snps", "K", "univ_p_endo", "univ_p_uf", "rho", "rho.lower",
         "rho.upper", "r2", "p", "p_fdr"]].sort_values(["CHR", "START"])
lv.columns = ["Block", "CHR", "Start (GRCh37)", "Stop (GRCh37)", "SNPs", "K", "Univariate P (endometriosis)",
              "Univariate P (fibroids)", "rho", "rho 95% CI lower", "rho 95% CI upper", "r2", "P", "FDR"]
panel("S1", "Panel D. Primary LAVA results with R9 case/control counts: all 88 valid bivariate tests", lv)

# ======================================================================= S2
m = pd.read_csv(RR / "MAGMA_cross_trait_comparison.csv")
# the comparison file is an inner join (18,300 genes); add genes tested in only one trait
TS = UF / "manuscript_assets_latest_2026_05_20/tables_supplementary"
a = pd.read_csv(TS / "S14a_MAGMA_Fibroids_AllGenes.tsv", sep="\t")
b = pd.read_csv(TS / "S14b_MAGMA_Endometriosis_AllGenes.tsv", sep="\t")
full = a.merge(b, on="GENE", how="outer", suffixes=("_f", "_e"))
for c in ("SYMBOL", "CHR", "START", "STOP"):
    full[c + "_f"] = full[c + "_f"].fillna(full[c + "_e"])
extra = full[~full.GENE.isin(m.GENE)].copy()
extra["group"] = np.where(extra.P_f < 2.5e-6, "significant only in fibroid analysis",
                 np.where(extra.P_e < 2.5e-6, "significant only in endometriosis analysis", "neither"))
extra.loc[extra.P_e.isna() | extra.P_f.isna(), "group"] = extra["group"] + " (tested in one trait only)"
m = pd.concat([m, extra], ignore_index=True)
cols = {"GENE": "Entrez ID", "SYMBOL_f": "Gene", "CHR_f": "CHR", "START_f": "Start (NCBI37)", "STOP_f": "Stop (NCBI37)",
        "NSNPS_f": "SNPs (fibroids)", "ZSTAT_f": "Z (fibroids)", "P_f": "P (fibroids)",
        "NSNPS_e": "SNPs (endometriosis)", "ZSTAT_e": "Z (endometriosis)", "P_e": "P (endometriosis)",
        "dz": "Z difference /sqrt(2)", "p_diff": "P (Z difference, approximate)", "group": "Significance group"}
mg = m[list(cols)].rename(columns=cols)
for tr, f in (("fibroids", "fibroids_neff"), ("endometriosis", "endo_neff")):
    p = RR / f"magma/{f}.genes.out"
    if p.exists():
        e = pd.read_csv(p, sep=r"\s+")[["GENE", "ZSTAT", "P"]].rename(
            columns={"GENE": "Entrez ID", "ZSTAT": f"Z (effective N, {tr})", "P": f"P (effective N, {tr})"})
        mg = mg.merge(e, on="Entrez ID", how="left")
mg = mg.sort_values("P (fibroids)")
panel("S2", "Panel A. MAGMA gene-level results (SNP-wise mean, NCBI Build 37.3, 1000 Genomes EUR). Bonferroni P < 2.5e-6. Groups are defined by threshold crossing and do not establish disorder specificity; the Z-difference test is approximate because Z scales with sample size.", mg)

loci = pd.read_csv(RR / "lava/lava_loci_2495_blocks.tsv", sep="\t")
valid = lv.set_index("Block")
rows = []
for _, g in m[m.group == "significant in both"].iterrows():
    b = loci[(loci.CHR == g.CHR_f) & (loci.STOP >= g.START_f) & (loci.START <= g.STOP_f)]
    for blk in b.LOC:
        if blk in valid.index:
            rows.append([g.SYMBOL_f, int(g.CHR_f), int(g.START_f), int(g.STOP_f), blk, valid.loc[blk, "rho"], valid.loc[blk, "FDR"],
                         "yes" if valid.loc[blk, "FDR"] < 0.05 else "no"])
        else:
            rows.append([g.SYMBOL_f, int(g.CHR_f), int(g.START_f), int(g.STOP_f), blk, np.nan, np.nan, "no valid test"])
panel("S2", "Panel B. The ten genes significant in both MAGMA analyses and the LAVA blocks overlapping them (primary analysis, R9 counts). Overlap reflects the same GWAS signal and is prioritisation, not replication.",
      pd.DataFrame(rows, columns=["Gene", "CHR", "Gene start (NCBI37.3)", "Gene stop (NCBI37.3)", "LAVA block", "rho", "FDR", "FDR-significant block"]))

# ======================================================================= S3
fm = pd.read_csv(RR / "finngen_finemap/FinnGen_insample_vs_our_SuSiE_by_locus.csv")
# 1000 Genomes-LD columns come from the repeat run at L = 10 (Panels D and F), so Panels A, D, F and Supplementary Figure S2 agree
_rr = pd.read_csv(RR / "susie/SuSiE_L_sensitivity_rerun.csv")
_rr = _rr[_rr.L == 10].assign(trait=lambda x: x.trait.replace({"Uterine_Fibroids": "Uterine fibroids"}))
_rm = pd.read_csv(RR / "susie/SuSiE_credible_set_membership_L10.csv").assign(trait=lambda x: x.trait.replace({"Uterine_Fibroids": "Uterine fibroids"}))
_rlead = _rm.sort_values("PIP", ascending=False).groupby(["locus", "trait", "credible_set"]).head(1)
_fmem = pd.read_csv(RR / "finngen_finemap/FinnGen_insample_credible_set_members_5loci.csv")
_inside = {k: int(g.BP.isin(set(_fmem[(_fmem.locus == k[0]) & (_fmem.trait == k[1])].position)).sum()) for k, g in _rlead.groupby(["locus", "trait"])}
fm["our_original_cs_L10"] = [int(_rr[(_rr.locus == l) & (_rr.trait == tr)].n_credible_sets.iloc[0]) for l, tr in zip(fm.locus, fm.trait)]
fm["converged_L10"] = ["yes" if bool(_rr[(_rr.locus == l) & (_rr.trait == tr)].converged.iloc[0]) else "no" for l, tr in zip(fm.locus, fm.trait)]
fm["our_leads_inside_finngen_cs"] = [_inside.get((l, tr), 0) for l, tr in zip(fm.locus, fm.trait)]
fm = fm[["locus", "trait", "finngen_cs_in_window", "finngen_good_cs", "finngen_cs_sizes", "finngen_lead_variants", "finngen_lead_gene",
         "finngen_lead_p", "our_original_cs_L10", "converged_L10", "our_leads_inside_finngen_cs"]]
fm["locus"] = fm["locus"].replace({"WT1": "chr11 (30.2 Mb)"})
panel("S3", "Panel A. Credible sets per locus: FinnGen R9 published SuSiE (in-sample LD; primary) versus SuSiE with the 1000 Genomes LD reference (L = 10; sensitivity analysis, same run as Panels D and F and Supplementary Figure S2). Only high-purity sets (good_cs) are counted in the text; sets were extracted when their lead variant lay within +/-500 kb of the regional GWAS lead variant. GRCh38.", fm.rename(columns={
    "locus": "Locus", "trait": "Trait", "our_original_cs_L10": "Credible sets, 1000G LD (L = 10)",
    "finngen_cs_in_window": "Credible sets, FinnGen in-sample", "finngen_good_cs": "High-purity sets (FinnGen good_cs)",
    "finngen_cs_sizes": "Set sizes", "finngen_lead_variants": "Lead variants", "finngen_lead_gene": "Gene of most severe consequence",
    "finngen_lead_p": "Lead P", "our_leads_inside_finngen_cs": "1000G-LD lead variants inside FinnGen sets", "converged_L10": "1000G LD run converged"}))
mem = pd.read_csv(RR / "finngen_finemap/FinnGen_insample_credible_set_members_5loci.csv")
mem["locus"] = mem["locus"].replace({"WT1": "chr11 (30.2 Mb)"})
mem = mem[["locus", "trait", "region", "cs", "v", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "p",
           "cs_specific_prob", "most_severe", "gene_most_severe"]]
panel("S3", "Panel B. Complete membership of FinnGen R9 in-sample credible sets in the five regions (GRCh38)", mem)
shared = []
for lc in mem.locus.unique():
    a = mem[(mem.locus == lc) & (mem.trait == "Uterine fibroids")]
    b = mem[(mem.locus == lc) & (mem.trait == "Endometriosis")]
    for (ra, ca), ga in a.groupby(["region", "cs"]):
        for (rb, cb), gb in b.groupby(["region", "cs"]):
            s = set(ga.v) & set(gb.v)
            if s:
                shared.append([lc, f"CS{ca} ({len(ga)} variants)", f"CS{cb} ({len(gb)} variants)", len(s), "; ".join(sorted(s)),
                               round(ga[ga.v.isin(s)].cs_specific_prob.sum(), 3), round(gb[gb.v.isin(s)].cs_specific_prob.sum(), 3)])
panel("S3", "Panel C. Credible sets shared between the two disorders (FinnGen in-sample). No shared variants at GREB1 or WNT4.",
      pd.DataFrame(shared, columns=["Locus", "Fibroid credible set", "Endometriosis credible set", "Shared variants",
                                    "Shared variant IDs", "Fibroid PIP in shared variants", "Endometriosis PIP in shared variants"]))
ls = pd.read_csv(RR / "susie/SuSiE_L_sensitivity_rerun.csv"); ls["locus"] = ls["locus"].replace({"WT1": "chr11 (30.2 Mb)"})
panel("S3", "Panel D. Sensitivity: 1000 Genomes-LD SuSiE (susie_rss, susieR 0.14.2; coverage 0.95, min_abs_corr 0.5; effective N) for L = 10, 5, 3, 1. Runs with converged = False did not reach convergence within 500 iterations.", ls)
dg = pd.read_csv(RR / "susie/SuSiE_LD_concordance_diagnostics.csv"); dg["locus"] = dg["locus"].replace({"WT1": "chr11 (30.2 Mb)"})
panel("S3", "Panel E. LD-concordance diagnostics for the 1000 Genomes-LD SuSiE: estimate_s_rss (0 = consistent) and number of variants flagged by kriging_rss (|z_std_diff| > 3 and logLR > 2)",
      dg.rename(columns={"estimated_s": "estimate_s_rss", "n_suspicious_variants_kriging": "Variants flagged (kriging_rss)"}))
cm = pd.read_csv(RR / "susie/SuSiE_credible_set_membership_L10.csv"); cm["locus"] = cm["locus"].replace({"WT1": "chr11 (30.2 Mb)"})
panel("S3", "Panel F. Sensitivity: credible-set membership of the 1000 Genomes-LD SuSiE, L = 10 (GRCh38 positions)", cm)

# ======================================================================= S4
old5 = pd.read_excel(OLD, sheet_name="Table S5 Single cell", header=None)
hdr = old5.index[old5[0] == "Gene"][0]
end = old5.index[old5[0].astype(str).str.startswith("Panel B")][0]
sc = old5.iloc[hdr + 1:end].copy()
sc.columns = old5.iloc[hdr].tolist()
sc = sc.dropna(how="all")
drop = ["Endometriosis_marker_celltype_enrichment_context", "Endometriosis_marker_celltype_enrichment_FDR",
        "Evidence_level", "Manuscript_safe_interpretation", "Guardrail"]
sc = sc.drop(columns=[c for c in drop if c in sc.columns])
sc.columns = ["Gene", "MAGMA -log10 P (fibroids)", "MAGMA -log10 P (endometriosis)",
              "Marker cell type, adenomyosis uterine tissue (GSE218044)", "Marker log2FC", "Marker FDR",
              "Fibroid tissue (GSE162122): top cell context", "Fibroid: mean log1p CP10K", "Fibroid: % cells expressing",
              "Myometrium (GSE162122): top cell context", "Myometrium: mean log1p CP10K", "Myometrium: % cells expressing"]
panel("S4", "Panel A. Single-cell context of the ten genes significant in both MAGMA analyses. Markers come from GSE218044 (uterine tissue from four individuals with adenomyosis; cell types assigned by cluster-level SingleR annotation in scEndoExplorer); fibroid and myometrium contexts come from GSE162122 (cell contexts assigned from marker signatures). Context only; no enrichment evidence is implied.", sc)
sl = pd.read_csv(RR / "SLDSC_full_statistics_from_raw_results.csv")
sl = sl[["trait", "baseline", "annotation", "prop_SNPs", "prop_h2", "enrichment", "enrichment_SE", "enrichment_P",
         "enrichment_FDR", "coef", "coef_SE", "coef_Z", "coef_P_onesided", "coef_FDR"]]
sl["annotation"] = sl["annotation"].str.replace("_", " ").str.replace("T cell  CD8 ", "T cell (CD8)", regex=False).str.strip()
sl.columns = ["Trait", "Baseline model", "Annotation", "Proportion of SNPs", "Proportion of h2", "Enrichment",
              "Enrichment SE", "Enrichment P", "Enrichment FDR", "Coefficient (tau)", "Coefficient SE",
              "Coefficient Z", "Coefficient P (one-sided)", "Coefficient FDR"]
panel("S4", "Panel B. S-LDSC for seven marker annotations under both baseline models. FDR is Benjamini-Hochberg over the seven annotations within trait and model. Figure 4 reports enrichment from baseline-LD v2.2 and the coefficient test from baseline v1.2.", sl)
panel("S4", "Panel C. Annotation construction", pd.DataFrame([
    ["Marker genes", "Seurat cluster markers with average log2 fold change > 0.25 and adjusted P < 0.05, ranked by adjusted P then fold change; top 500 per cell type (B cell 419, CD8 T cell 496)"],
    ["Source", "GSE218044: uterine tissue from four individuals with adenomyosis (Chen et al. 2022), 39,995 cells after subsampling, processed in scEndoExplorer (cell_type_marker_genes.csv)"],
    ["Cell labels", "Cluster-level SingleR annotation (Human Primary Cell Atlas and Blueprint references), harmonised to seven uterine cell classes. A cell-level reannotation with confidence scores is reported in Panels E-H."],
    ["Gene window", "Gene body +/- 100 kb"],
    ["Reference panel", "1000 Genomes Phase 3 European; HapMap3 regression SNPs"],
    ["Baseline models", "baseline v1.2 (coefficient test) and baseline-LD v2.2 (enrichment)"],
    ["Software", "LDSC v1.0.1 (Python 3 fork)"],
], columns=["Item", "Specification"]))
panel("S4", "Panel D. Overlap between marker gene sets", pd.read_csv(RR / "SLDSC_marker_gene_set_overlap.csv"))
RA = RR / "reannotation"
ct = pd.read_csv(RA / "reannotation_summary_by_celltype.csv")
panel("S4", "Panel E. Cell-level reannotation of GSE218044 (SingleR per cell, HPCA and Blueprint): agreement with cluster-level labels by cell class. Concordant = per-cell label equals cluster label and not flagged; pruned = SingleR outlier flag; delta = score margin to the next-best label.",
      ct.rename(columns={"cluster_label": "Cluster-level label", "cells": "Cells", "concordant": "Fraction concordant", "discordant": "Fraction discordant", "pruned": "Fraction pruned", "median_delta": "Median delta.next"}))
dn = pd.read_csv(RA / "reannotation_summary_by_donor.csv").sort_values(["cluster_label", "donor"])
panel("S4", "Panel F. Cell-level reannotation by donor and cell class", dn.rename(columns={"donor": "Donor", "cluster_label": "Cluster-level label", "cells": "Cells", "concordant": "Fraction concordant", "discordant": "Fraction discordant", "pruned": "Fraction pruned"}))
cf = pd.read_csv(RA / "reannotation_confusion_cluster_vs_cell.csv")
panel("S4", "Panel G. Cluster-level label (rows) versus per-cell label (columns), number of cells", cf.rename(columns={"cluster_label": "Cluster-level label"}))
ov2 = pd.read_csv(RA / "marker_overlap_original_vs_confident.csv")
panel("S4", "Panel H(i). Marker genes re-derived from concordantly labelled cells (23,587 cells) versus the original markers. The procedure reproduced the original marker lists exactly when applied to all cells. Re-derived markers for smooth muscle and B cells rest on few concordant cells (231 and 79), so these two sets are less reliable.", ov2)
sg = pd.read_csv(RA / "shared_gene_marker_status.csv").fillna("")
panel("S4", "Panel H(ii). Marker status of the ten shared MAGMA genes, original versus re-derived markers", sg)
cmp = pd.read_csv(RR / "sldsc_reannot/SLDSC_original_vs_confident_markers.csv")
cmp["ct"] = cmp.ct.str.replace("T_cell__CD8_", "CD8 T cell", regex=False).str.replace("_", " ")
cmp["model"] = cmp.model.map({"v12": "baseline v1.2", "ld22": "baseline-LD v2.2"})
cmp["set"] = cmp.set.map({"original": "original markers", "confident": "re-derived markers (concordant cells)"})
panel("S4", "Panel H(iii). S-LDSC with original versus re-derived markers, both baseline models. FDR over seven annotations within trait, model and marker set.",
      cmp.rename(columns={"set": "Marker set", "trait": "Trait", "model": "Baseline model", "ct": "Annotation", "prop": "Proportion of SNPs", "enr": "Enrichment", "enr_p": "Enrichment P", "enr_fdr": "Enrichment FDR", "z": "Coefficient Z", "coef_p": "Coefficient P (one-sided)", "coef_fdr": "Coefficient FDR"})[
          ["Marker set", "Trait", "Baseline model", "Annotation", "Proportion of SNPs", "Enrichment", "Enrichment P", "Enrichment FDR", "Coefficient Z", "Coefficient P (one-sided)", "Coefficient FDR"]])

import glob as _glob, os as _os
from scipy.stats import norm as _norm
_rows = []
for _set in ("original", "confident"):
    for _tr in ("fibroids", "endometriosis"):
        for _m in ("v12", "ld22"):
            _d = RR / f"sldsc_reannot/results_control_{_set}_{_tr}_{_m}"
            _t = []
            for _f in sorted(_glob.glob(str(_d / "*.results"))):
                _r = pd.read_csv(_f, sep="\t").iloc[-1]
                _t.append([_os.path.basename(_f)[:-8], _r["Prop._SNPs"], _r.Enrichment, _r.Enrichment_p, _r["Coefficient_z-score"]])
            if not _t:
                continue
            _t = pd.DataFrame(_t, columns=["Annotation", "Proportion of SNPs", "Enrichment", "Enrichment P", "Coefficient Z"])
            _t["Coefficient P (one-sided)"] = _norm.sf(_t["Coefficient Z"])
            _t["Enrichment FDR"] = multipletests(_t["Enrichment P"], method="fdr_bh")[1]
            _t["Coefficient FDR"] = multipletests(_t["Coefficient P (one-sided)"], method="fdr_bh")[1]
            _t.insert(0, "Baseline model", {"v12": "baseline v1.2", "ld22": "baseline-LD v2.2"}[_m])
            _t.insert(0, "Trait", _tr)
            _t.insert(0, "Marker set", {"original": "original markers", "confident": "re-derived markers (concordant cells)"}[_set])
            _rows.append(_t)
if _rows:
    _ctl = pd.concat(_rows, ignore_index=True)
    _ctl["Annotation"] = _ctl["Annotation"].str.replace("T_cell__CD8_", "CD8 T cell", regex=False).str.replace("_", " ")
    panel("S4", "Panel I. S-LDSC with an all-genes control annotation (all GSE218044 genes with coordinates, +/-100 kb) fitted jointly with the baseline model and each cell-type annotation. FDR over seven annotations within trait, model and marker set.", _ctl)

# ======================================================================= S5
o = pd.read_csv(UF / "outputs_v2/formal_ld_audit/lava/Table_LAVA_Full_2495Blocks_AllResults_overlap_ldsc_standardized_withFDR.csv")
orig16 = list(o[o.p_fdr < 0.05].sort_values("p").LOC)
runs = {"0 (no overlap)": "res_grid_corr_0.0.csv", "0.1724 (primary)": "LAVA_corrected_N_ov01724_all_bivariate_rows.csv",
        "0.259 (overlap count, munged-file totals)": "res_grid_corr_0.2590.csv", "0.35": "res_grid_corr_0.35.csv",
        "0.369 (overlap count, Release 9 totals)": "res_grid_corr_0.3685.csv"}
summ, wide = [], o[o.LOC.isin(orig16)][["LOC", "CHR", "START", "STOP", "rho", "p_fdr"]].rename(
    columns={"LOC": "Block", "START": "Start (GRCh37)", "STOP": "Stop (GRCh37)", "rho": "rho (munged-file counts)", "p_fdr": "FDR (munged-file counts)"})
allv = {}
for lab, f in runs.items():
    d = pd.read_csv(RR / "lava" / f)
    v = d[d.p.notna()].copy()
    v["fdr"] = multipletests(v.p, method="fdr_bh")[1]
    allv[lab] = (d, v)
    x = d[d.LOC.isin(orig16)].merge(v[["LOC", "fdr"]], on="LOC", how="left")
    summ.append([lab, len(d), len(v), int((v.fdr < 0.05).sum()), int(x.p.notna().sum()), int((x.fdr < 0.05).sum()),
                 int((x.rho > 0).sum()), ", ".join(sorted(set(v[v.fdr < 0.05].LOC) - set(orig16))) or "none"])
    x = x[["LOC", "rho", "rho.lower", "rho.upper", "fdr"]].rename(columns={
        "LOC": "Block", "rho": f"rho [{lab}]", "rho.lower": f"CI lower [{lab}]", "rho.upper": f"CI upper [{lab}]", "fdr": f"FDR [{lab}]"})
    wide = wide.merge(x, on="Block", how="left")
rep = pd.read_csv(RR / "lava/res_repro_legacy.csv")
rv = rep[rep.p.notna()].copy(); rv["fdr"] = multipletests(rv.p, method="fdr_bh")[1]
summ.append(["0.1724, case/control counts of the munged input files", len(rep), len(rv), int((rv.fdr < 0.05).sum()), 16,
             int((rv[rv.LOC.isin(orig16)].fdr < 0.05).sum()), 16, "none"])
panel("S5", "Panel A. LAVA results by sampling correlation, R9 case/control counts, same 150 candidate blocks. The final row uses the case/control counts carried in the munged input files (35,474/183,254 and 17,045/222,290).",
      pd.DataFrame(summ, columns=["Sampling correlation", "Bivariate rows", "Valid tests", "FDR-significant blocks",
                                  "Primary 16: estimable", "Primary 16: FDR-significant", "Primary 16: positive rho",
                                  "Significant blocks outside the primary 16"]))
stable = {"block_1074", "block_1075", "block_1076", "block_1077", "block_0024", "block_0216", "block_0666", "block_1644"}
wide.insert(1, "Classification", wide.Block.map(lambda b: "stable core" if b in stable else "sensitivity-dependent"))
wide.insert(2, "Region", wide.Block.map({"block_1074": "chr6 ESR1/SYNE1", "block_1075": "chr6 ESR1/SYNE1", "block_1076": "chr6 ESR1/SYNE1",
                                          "block_1077": "chr6 ESR1/SYNE1", "block_0024": "chr1 WNT4", "block_0216": "chr2 GREB1",
                                          "block_0666": "chr4 70.5 Mb", "block_1642": "chr11 29.5 Mb", "block_1644": "chr11 32.1 Mb"}).fillna(""))
panel("S5", "Panel B. The 16 primary blocks across sampling correlations: rho, 95% CI and FDR (FDR within each run). Blank = not estimable. Stable core = FDR < 0.05 at 0.1724, 0.259, 0.35 and 0.369, and whenever estimable at 0; adjacent blocks form one region.",
      wide.sort_values(["Classification", "CHR", "Start (GRCh37)"], ascending=[False, True, True]))

d0, v0 = allv["0 (no overlap)"]
_, vp = allv["0.1724 (primary)"]
dp = allv["0.1724 (primary)"][0]
extra = v0[(v0.fdr < 0.05) & ~v0.LOC.isin(orig16)][["LOC", "CHR", "START", "STOP", "rho", "rho.lower", "rho.upper", "p", "fdr"]]
extra = extra.merge(dp[["LOC", "rho", "p"]].rename(columns={"rho": "rho [0.1724 (primary)]", "p": "P [0.1724 (primary)]"}), on="LOC", how="left")
extra = extra.merge(vp[["LOC", "fdr"]].rename(columns={"fdr": "FDR [0.1724 (primary)]"}), on="LOC", how="left")
panel("S5", "Panel C. Blocks significant only without overlap modelling (sampling correlation 0): rho, 95% CI, P and FDR in that run, with the estimate from the primary analysis for comparison. A rho of 1 is the parameter bound.",
      extra.rename(columns={"LOC": "Block", "START": "Start (GRCh37)", "STOP": "Stop (GRCh37)", "rho": "rho [0 (no overlap)]",
                            "rho.lower": "CI lower [0 (no overlap)]", "rho.upper": "CI upper [0 (no overlap)]",
                            "p": "P [0 (no overlap)]", "fdr": "FDR [0 (no overlap)]"}).sort_values(["CHR", "Start (GRCh37)"]))

# ======================================================================= S6
panel("S6", "Panel A. Genome build of each input", pd.DataFrame([
    ["FinnGen R9 GWAS summary statistics", "GRCh38"], ["FinnGen R9 SuSiE fine-mapping", "GRCh38"],
    ["SuSiE windows and lead positions (Table S3)", "GRCh38"], ["1000 Genomes Phase 3 EUR PLINK reference (LAVA, SuSiE LD)", "GRCh37"],
    ["LAVA block definitions (2,495 blocks)", "GRCh37"], ["MAGMA gene boundaries (NCBI 37.3)", "GRCh37"],
    ["LDSC reference LD scores", "rsID-based (HapMap3)"],
    ["Harmonisation", "Variants matched across builds by rsID and alleles; no liftOver of summary statistics"],
], columns=["Component", "Build"]))
raw = pd.read_csv(RR / "build_map/lead_variants_raw.csv")
rs = pd.read_csv(RR / "build_map/finngen_rsids.csv").drop_duplicates(["chr", "pos38", "rsids"])
bim = pd.read_csv(RR / "build_map/bim_hits.txt", sep="\t", header=None, names=["c", "rsid", "cm", "pos37", "a1", "a2"])
raw["chr"] = raw.chr.astype(str); rs["chr"] = rs.chr.astype(str)
mp = raw.merge(rs.groupby(["chr", "pos38"]).agg(ref=("ref", "first"), alt=("alt", lambda x: "/".join(x)), rsid=("rsids", "first")).reset_index(),
               on=["chr", "pos38"], how="left").merge(bim[["rsid", "pos37"]], on="rsid", how="left")
def block(r):
    if pd.isna(r.pos37): return "not in 1000G reference"
    b = loci[(loci.CHR == int(r.chr)) & (loci.START <= r.pos37) & (loci.STOP >= r.pos37)]
    return b.LOC.iloc[0] if len(b) else ""
mp["LAVA block (GRCh37)"] = mp.apply(block, axis=1)
mp["Block FDR (primary)"] = mp["LAVA block (GRCh37)"].map(lambda b: valid.loc[b, "FDR"] if b in valid.index else np.nan)
mp = mp[mp.locus.isin(["ESR1_SYNE1", "WNT4", "GREB1", "WT1"])]
mp["locus"] = mp.locus.replace({"WT1": "chr11 (30.2 Mb)"})
mp = mp[["locus", "trait", "cs_id", "rsid", "chr", "ref", "alt", "pos38", "pos37", "LAVA block (GRCh37)", "Block FDR (primary)", "p", "cs_size"]]
mp.columns = ["Locus", "Source", "FinnGen credible set", "rsID", "CHR", "Ref", "Alt", "Position GRCh38", "Position GRCh37",
              "LAVA block (GRCh37)", "Block FDR (primary)", "GWAS P", "Credible-set size"]
panel("S6", "Panel B. Cross-build mapping of FinnGen credible-set lead variants and SuSiE window leads at ESR1/SYNE1, WNT4, GREB1 and chr11 (30.2 Mb). Indels absent from the 1000 Genomes reference by rsID are marked.", mp)
panel("S6", "Panel C. MAGMA gene boundaries (NCBI 37.3) of the ten shared genes and the LAVA blocks they overlap: seven lie in FDR-significant blocks",
      sheets["S2"][1][1])

# ======================================================================= write
titles = {"S1": "Data sources, sample sizes, genome-wide and primary local genetic correlation",
          "S2": "MAGMA gene-level association", "S3": "Fine-mapping", "S4": "S-LDSC and single-cell context",
          "S5": "LAVA sensitivity to the sampling correlation", "S6": "Genome-build traceability"}
names = {k: f"Table {k} " + {"S1": "Sample size LDSC LAVA", "S2": "MAGMA", "S3": "Fine mapping", "S4": "S-LDSC single cell",
                              "S5": "LAVA sensitivity", "S6": "Genome builds"}[k] for k in titles}
with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    idx = pd.DataFrame([[f"Supplementary Table {k}", names[k], titles[k]] for k in titles],
                       columns=["Supplementary table", "Worksheet", "Title"])
    pd.DataFrame([["Supplementary Tables: Integrative genomic analyses of uterine fibroids and endometriosis"]]).to_excel(
        xw, sheet_name="Index", header=False, index=False)
    idx.to_excel(xw, sheet_name="Index", startrow=2, index=False)
    for k, panels in sheets.items():
        r = 0
        for title, df in panels:
            pd.DataFrame([[title]]).to_excel(xw, sheet_name=names[k], startrow=r, header=False, index=False)
            df.to_excel(xw, sheet_name=names[k], startrow=r + 1, index=False)
            r += len(df) + 4
import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill
from openpyxl.utils import get_column_letter

# Plain formatting matching the originally submitted workbook: Calibri 11, no fills,
# no borders, no frozen panes; only panel titles and column-header rows are bold.
PLAIN, BOLD = Font(name="Calibri", size=11), Font(name="Calibri", size=11, bold=True)
wb = openpyxl.load_workbook(OUT)
for ws in wb.worksheets:
    ws.freeze_panes = None
    bold_rows = set()
    for row in ws.iter_rows(min_col=1, max_col=1):
        c = row[0]
        if isinstance(c.value, str) and (c.value.startswith("Panel") or c.value.startswith("Supplementary Tables:")):
            bold_rows.update({c.row, c.row + 1})
    if ws.title == "Index":
        bold_rows.update({1, 3})
    widths = {}
    for row in ws.iter_rows():
        for c in row:
            c.font = BOLD if c.row in bold_rows else PLAIN
            c.border = Border()
            c.fill = PatternFill(fill_type=None)
            c.alignment = Alignment()
            if isinstance(c.value, float) and c.value != 0 and abs(c.value) < 1e-3:
                c.number_format = "0.00E+00"
            is_title = isinstance(c.value, str) and c.value.startswith(("Panel", "Supplementary Tables:"))
            if c.value is not None and not is_title and c.row < 300:
                widths[c.column] = max(widths.get(c.column, 0), len(str(c.value)))
    for col, w in widths.items():
        ws.column_dimensions[get_column_letter(col)].width = min(max(w + 2, 10), 45)
wb.save(OUT)
print("wrote", OUT)
for k, panels in sheets.items():
    print(k, [(t[:40], len(d)) for t, d in panels])
