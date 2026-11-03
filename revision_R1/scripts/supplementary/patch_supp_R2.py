s=open('build_supplementary_R1.py').read()
def rep(a,b,cnt=1):
    global s
    assert s.count(a)==cnt,(s.count(a),a[:60]); s=s.replace(a,b)
rep('Supplementary_Tables_S1_to_S6_F_IG_R1.xlsx','Supplementary_Tables_S1_to_S6_F_IG_R2.xlsx',2)
rep('MAGMA: munged-file N (same significant genes with R9 effective N, Table S2)','MAGMA: R9 effective N (primary; munged-file N as sensitivity, same significant genes, Table S2)',2)
rep('"LDSC and MAGMA (constant N in munged files)"','"LDSC and MAGMA, first run (constant N in munged files)"',2)
old=s[s.index('for tr, f in (("fibroids", "fibroids_neff"), ("endometriosis", "endo_neff")):'):s.index('mg = mg.sort_values("P (fibroids)")')]
new='''# R2: primary statistics = MAGMA with R9 effective N; the first run (munged-file N) is kept as sensitivity
mg = mg.rename(columns={"Z (fibroids)": "Z (fibroids, munged-file N)", "P (fibroids)": "P (fibroids, munged-file N)",
                        "Z (endometriosis)": "Z (endometriosis, munged-file N)", "P (endometriosis)": "P (endometriosis, munged-file N)"})
for tr, f, ab in (("fibroids", "fibroids_neff", "fibroids"), ("endometriosis", "endo_neff", "endometriosis")):
    e = pd.read_csv(RR / f"magma/{f}.genes.out", sep=r"\\s+")[["GENE", "ZSTAT", "P"]].rename(
        columns={"GENE": "Entrez ID", "ZSTAT": f"Z ({ab})", "P": f"P ({ab})"})
    mg = mg.merge(e, on="Entrez ID", how="left")
mg = mg[["Entrez ID", "Gene", "CHR", "Start (NCBI37)", "Stop (NCBI37)", "SNPs (fibroids)", "Z (fibroids)", "P (fibroids)",
         "SNPs (endometriosis)", "Z (endometriosis)", "P (endometriosis)", "Significance group",
         "Z (fibroids, munged-file N)", "P (fibroids, munged-file N)", "Z (endometriosis, munged-file N)", "P (endometriosis, munged-file N)"]]
'''
s=s.replace(old,new)
rep('mg = mg.sort_values("P (fibroids)")','mg = mg.sort_values("P (fibroids)", na_position="last")')
rep('panel("S2", "Panel A. MAGMA gene-level results (SNP-wise mean, NCBI Build 37.3, 1000 Genomes EUR). Bonferroni P < 2.5e-6.',
    'panel("S2", "Panel A. MAGMA gene-level results (SNP-wise mean, NCBI Build 37.3, 1000 Genomes EUR) with the Release 9 effective sample sizes (primary; files magma/fibroids_neff.genes.out and magma/endo_neff.genes.out; source of Figure 3). The last four columns give the first run with the constant sample sizes of the munged files (sensitivity; same significant genes). Bonferroni P < 2.5e-6.')
rep('"Panel B. S-LDSC for seven marker annotations under both baseline models. FDR is Benjamini-Hochberg over the seven annotations within trait and model. Figure 4 reports enrichment from baseline-LD v2.2 and the one-sided coefficient test (prespecified positive enrichment) from baseline v1.2; two-sided values are given for comparison."',
    '"Panel B. S-LDSC for seven marker annotations under both baseline models (original markers). FDR is Benjamini-Hochberg over the seven annotations within trait and model. Figure 4 reports enrichment from baseline-LD v2.2 and the one-sided coefficient test (positive enrichment, following Finucane et al.; there is no record that the direction was fixed before these results were examined) from baseline v1.2; two-sided values are given for comparison and in Panel H(iv)."')
rep('Only high-purity sets (good_cs) are counted in the text; sets were extracted when their lead variant lay within +/-500 kb of the regional GWAS lead variant. GRCh38.',
    'Only high-purity sets (good_cs) are counted in the text; in the downloaded files good_cs is true exactly when the minimum pairwise r2 within the set is at least 0.25 (|r| >= 0.5). Sets were extracted when their lead variant lay within +/-500 kb of the regional GWAS lead variant. Source files: finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.SUSIE.cred.summary.tsv and finngen_R9_N14_ENDOMETRIOSIS.SUSIE.cred.summary.tsv (with the matching .snp.filter.tsv files). GRCh38.')
anchor='import glob as _glob, os as _os'
add='''dd = pd.read_csv(RR / "sldsc_R2/SLDSC_direction_and_marker_dependence.csv")
dd["ct"] = dd.ct.str.replace("T_cell__CD8_", "CD8 T cell", regex=False).str.replace("_", " ")
dd["model"] = dd.model.map({"v12": "baseline v1.2", "ld22": "baseline-LD v2.2"})
dd = dd.rename(columns={"set": "Marker set", "trait": "Trait", "model": "Baseline model", "ct": "Annotation", "enr": "Enrichment", "enr_fdr": "Enrichment FDR",
                        "z": "Coefficient Z", "coef_p": "Coefficient P (one-sided)", "coef_p_two": "Coefficient P (two-sided)",
                        "coef_fdr_one_sided": "Coefficient FDR (one-sided)", "coef_fdr_two_sided": "Coefficient FDR (two-sided)"})
panel("S4", "Panel H(iv). Dependence of the conditional coefficient result on marker assignment (original versus re-derived markers) and test direction (one-sided versus two-sided). Two-sided P = 2 x min(p, 1 - p) of the one-sided P. FDR is Benjamini-Hochberg over the seven annotations within trait, baseline model and marker set. The one-sided direction (positive enrichment) follows the framework of Finucane et al. and there is no record that it was fixed before these results were examined.",
      dd[["Marker set", "Trait", "Baseline model", "Annotation", "Enrichment", "Enrichment FDR", "Coefficient Z", "Coefficient P (one-sided)",
          "Coefficient P (two-sided)", "Coefficient FDR (one-sided)", "Coefficient FDR (two-sided)"]])
_RA = Path("/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/github_repo/revision_R1/results/reannotation")
_mk = []
for _f, _lab in (("markers_all_cells_reproduction.csv", "original markers (all cells; reproduced)"), ("markers_confident_cells.csv", "re-derived markers (concordant cells)")):
    _x = pd.read_csv(_RA / _f)
    _x.insert(0, "Marker set", _lab)
    _x["rank"] = _x.groupby("cell_type").cumcount() + 1
    _mk.append(_x)
panel("S4", "Panel J. Marker gene lists used for the S-LDSC annotations (an adjusted P of 0 is numerical underflow, below about 1e-308, not an exact zero): top 500 markers per cell class (fewer where fewer passed the thresholds), ranked by adjusted P and then fold change. Cell-level reannotation code: scripts/reannotation/reannotate_GSE218044.R in the code repository; per-cell labels: results/reannotation/GSE218044_cell_level_labels.csv.gz.",
      pd.concat(_mk, ignore_index=True)[["Marker set", "cell_type", "rank", "gene", "avg_log2FC", "p_val_adj"]].rename(
          columns={"cell_type": "Cell class", "rank": "Rank", "gene": "Gene", "avg_log2FC": "Average log2 fold change", "p_val_adj": "Adjusted P"}))

'''
rep(anchor,add+anchor)
anchor2='# ======================================================================= S6'
add2='''LR2 = RR / "lava_R2"
lg = pd.read_csv(LR2 / "LAVA_all148_rows_by_specification_long.csv")
lg["START"] = lg.START.astype("Int64"); lg["STOP"] = lg.STOP.astype("Int64")
lg = lg.rename(columns={"specification": "Sampling correlation", "LOC": "Block", "START": "Start (GRCh37)", "STOP": "Stop (GRCh37)", "n_snps": "SNPs",
                        "univ_p_endo": "Univariate P (endometriosis)", "univ_p_uf": "Univariate P (fibroids)", "h2_endo": "Local h2 (endometriosis)",
                        "h2_uf": "Local h2 (fibroids)", "rho.lower": "rho 95% CI lower", "rho.upper": "rho 95% CI upper", "p": "P", "valid_P": "Valid P",
                        "fdr_within_run": "FDR within run", "fdr_fixed88": "FDR, fixed family (BH, m = 88)", "fdr_fixed148": "FDR, fixed family (BH, m = 148; no valid P = 1)", "reason_no_valid_P": "Reason no valid P"})
panel("S5", "Panel D. All 148 bivariate rows under each of the five sampling-correlation specifications (740 rows): raw P, estimate, 95% CI, FDR within the run, and FDR in a fixed correction family (Benjamini-Hochberg with m = 88, the number of valid tests in the primary analysis; only blocks with a valid P in the primary analysis enter this column), and in the most conservative fixed family (m = 148, all candidate rows enter and rows without a valid P are assigned P = 1). LAVA returns the bivariate estimate as NA when it is outside the bound (|rho| > 1.25) or not estimable and does not report which; such rows carry no P value. The two candidate blocks that failed locus processing in every run are listed in Panel E.",
      lg)
fl = pd.read_csv(LR2 / "LAVA_failed_locus_processing.csv").rename(columns={"LOC": "Block", "status": "Status", "msg": "Message"})
panel("S5", "Panel E. The two candidate blocks that failed during locus processing in every specification", fl)
st = pd.read_csv(LR2 / "LAVA_stable_core_three_definitions.csv")
keep = ["Block", "CHR", "START", "STOP"]
for k in ("0", "0.1724", "0.259", "0.35", "0.369"):
    keep += [f"estimable [{k}]", f"FDR within run [{k}]", f"FDR fixed family m=88 [{k}]", f"FDR fixed family m=148 [{k}]", f"FDR common set [{k}]"]
keep += ["stable: within-run FDR", "stable: fixed family", "stable: fixed family m=148", "tested without overlap", "in common estimable set"]
st = st[keep].rename(columns={"START": "Start (GRCh37)", "STOP": "Stop (GRCh37)"})
panel("S5", "Panel F. The 16 primary blocks classified under four correction schemes: FDR within each run (as in the main analysis), a fixed family of the 88 primary tests, the most conservative fixed family of all 148 candidate rows (rows without a valid P assigned P = 1) and the 14 blocks estimable in all five runs (BH within that set). The stable core is eight blocks in five regions with FDR within each run, seven blocks in four regions with the 88-test family (chromosome 11 at 32.1 Mb drops out) and six blocks in four regions with the 148-row family (the fourth ESR1/SYNE1 block, block_1077, also drops out; its FDR is 0.051 at 0.369). Stable = FDR < 0.05 at 0.1724, 0.259, 0.35 and 0.369 and, where estimable, at 0. Blocks not estimable without overlap (block_0216, block_0666, block_1077) were not tested successfully under that assumption; their stability refers to the overlap-aware values only. Eight of the 16 blocks lie outside the common estimable set.", st)
_info = pd.DataFrame([
    ["Intercept, fibroids (LDSC, effective N)", 1.0967, "diagonal element of the intercept matrix"],
    ["Intercept, endometriosis (LDSC, effective N)", 1.1274, "diagonal element"],
    ["Cross-trait intercept (LDSC)", 0.1917, "off-diagonal element (SE about 0.009)"],
    ["Standardised sampling correlation = 0.1917 / sqrt(1.0967 x 1.1274)", 0.1724, "primary value; approximate delta-method SE about 0.008 from the intercept SEs (95% interval about 0.157 to 0.188), ignoring covariances between the three intercept estimates"],
    ["Count-derived, munged-file totals: 59,264 / sqrt(218,728 x 239,335)", 0.259, "depends on an overlap count that could not be verified"],
    ["Count-derived, Release 9 totals: 59,264 / sqrt(210,870 x 122,652)", 0.369, "same unverified count"],
    ["Intermediate value", 0.35, "bracketing value"],
], columns=["Quantity", "Value", "Note"])
panel("S5", "Panel G. Derivation of the sampling correlations. The LDSC cross-trait intercept reflects the covariance of the Z-score errors of the two GWAS induced by shared participants (shared cases and shared controls raise it; a participant who is a case in one study and a control in the other lowers it); dividing by the square root of the product of the two intercepts gives the correlation of the sampling errors that LAVA requires. The intercept also absorbs any shared confounding, so 0.1724 is a model-based approximation of the overlap, not a count.", _info)

'''
rep(anchor2,add2+anchor2)
rep('"Panel B. The 16 primary blocks across sampling correlations: rho, 95% CI and FDR (FDR within each run). Blank = not estimable. Stable core = FDR < 0.05 at 0.1724, 0.259, 0.35 and 0.369, and whenever estimable at 0; adjacent blocks form one region."',
    '"Panel B. The 16 primary blocks across sampling correlations: rho, 95% CI and FDR (FDR within each run). Blank = not estimable (no valid P in that run). Stable core = FDR < 0.05 at 0.1724, 0.259, 0.35 and 0.369, and whenever estimable at 0; three of the eight stable blocks (block_0216, block_0666, block_1077) have no estimate without overlap. Adjacent blocks form one region. Alternative correction families are in Panel F."')
rep('["Estimate outside parameter bound (|rho| > 1.25), P set to missing", 58, 60],','["No valid P (estimate returned as NA: outside |rho| <= 1.25 or not estimable; LAVA does not report which)", 58, 60],')
rep('f"CS{ca} ({len(ga)} variants)", f"CS{cb} ({len(gb)} variants)"','f"CS{ca} ({len(ga)} variant{\'s\' if len(ga) != 1 else \'\'})", f"CS{cb} ({len(gb)} variant{\'s\' if len(gb) != 1 else \'\'})"')
rep('["Intercept, fibroids (LDSC, effective N)", 1.0967, "diagonal element of the intercept matrix"]','["Intercept, fibroids (GenomicSEM LDSC, munged-file sample sizes of the earlier configuration)", 1.0967, "diagonal element of the intercept matrix"]')
rep('["Intercept, endometriosis (LDSC, effective N)", 1.1274, "diagonal element"]','["Intercept, endometriosis (GenomicSEM LDSC, munged-file sample sizes of the earlier configuration)", 1.1274, "diagonal element"]')
rep('["Cross-trait intercept (LDSC)", 0.1917, "off-diagonal element (SE about 0.009)"]','["Cross-trait intercept (GenomicSEM LDSC, munged-file sample sizes of the earlier configuration)", 0.1917, "off-diagonal element (SE about 0.009); the fit with the Release 9 effective sample sizes has cross-trait intercept 0.1906 and gives a standardised value of 0.1716"]')
rep("(with the matching .snp.filter.tsv files). GRCh38.","(with the matching .snp.filter.tsv files); in the 141 sets in these files the smallest flagged r2 is 0.251 and the largest unflagged r2 is 0.243. GRCh38.")

a3 = '# ======================================================================= write'
assert s.count(a3) == 1
s = s.replace(a3, '_j = [p for p in sheets["S4"] if p[0].startswith("Panel J")]\nsheets["S4"] = [p for p in sheets["S4"] if not p[0].startswith("Panel J")] + _j\n' + a3)
open('build_supplementary_R2.py','w').write(s)
