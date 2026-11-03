s=open('edit_manuscript_R2.py').read()
def r(a,b,cnt=1):
    global s
    assert s.count(a)==cnt,(s.count(a),a[:70]); s=s.replace(a,b)
r(" Two-sided FDRs can be smaller than one-sided ones, because the two-sided family also counts annotations with strongly negative coefficients as rejections (with the re-derived markers, the endometriosis stromal coefficient FDR was 0.041 two-sided and 0.062 one-sided; the fibroid values were 0.21 and 0.11).","")
r("(stromal coefficient FDR with the original markers, two-sided: 0.090 for fibroids","(stromal coefficient FDR with the original markers, two-sided: 0.089 for fibroids")
r('"eight of which, across five genomic regions, remained significant under all overlap-aware sampling correlations")','"eight of which, across five genomic regions, remained significant under all overlap-aware sampling correlations when FDR was calculated within each run (six or seven blocks in four regions under fixed correction families)")')
r('"were robust to assumptions about sample overlap", "remained significant under all overlap-aware sampling correlations")','"were robust to assumptions about sample overlap", "remained significant under all overlap-aware sampling correlations with FDR calculated within each run")')
r("all 148 bivariate rows in each specification with raw P values and reasons for missing values,","all 148 bivariate rows in each specification with raw P values and missing-result flags,")
r("raw P values, estimates and reasons for missing P values for all 148 rows in each specification","raw P values, estimates and missing-result flags for all 148 rows in each specification")
r('NEW = {"3":','NEW = {"2": FIG + "revision_figures/Figure_2_LAVA_R2.png", "3":')
extra='''
# ---------------------------------------------------------------- audit fixes
edit("Cell-type heritability enrichment analysis using S-LDSC", "endometriosis z = 2.53, FDR = 0.041; baseline v1.2", "endometriosis z = 2.53, FDR = 0.040; baseline v1.2")
edit("Cell-type heritability enrichment analysis using S-LDSC", "would be 0.090 for fibroids", "would be 0.089 for fibroids")
edit("Cell-type heritability enrichment analysis using S-LDSC", "endometriosis z = 2.37, FDR = 0.062).",
     "endometriosis z = 2.37, FDR = 0.062). Two-sided FDRs can be smaller than one-sided ones, because the two-sided family also counts annotations with strongly negative coefficients as rejections (with the re-derived markers, the endometriosis stromal coefficient FDR was 0.041 two-sided and 0.062 one-sided; the fibroid values were 0.21 and 0.11).")
edit("Cell-type heritability enrichment analysis using S-LDSC", "may not fully capture. Likewise, the absence of smooth-muscle enrichment does not contradict the smooth-muscle histology of fibroids [6, 16].",
     "may not fully capture [6, 16]. Likewise, the absence of smooth-muscle enrichment does not contradict the smooth-muscle histology of fibroids [2].")
edit("Local genetic correlation was estimated with LAVA", "and estimates outside the parameter bound (|ρ| > 1.25) were set to missing.",
     "and no estimate was returned when it lay outside the parameter bound (|ρ| > 1.25) or could not be estimated.")
edit("Local genetic correlation was estimated with LAVA", "was derived from the LDSC intercept covariance matrix: diagonal",
     "was derived from the LDSC intercept covariance matrix of the fit with the munged-file sample sizes of the earlier configuration: diagonal")
edit("The values 0.259 and 0.369 were calculated", "both values were used only in sensitivity analyses.",
     "both values were used only in sensitivity analyses. The value 0.35 was included as an intermediate value between them.")
edit("Publicly available single-cell RNA-seq data were used", "so two-sided values are also reported (Supplementary Table S4).",
     "so two-sided values are also reported (Supplementary Table S4). Coefficient tests were also computed under baseline-LD v2.2 (Supplementary Table S4, Panel B).")
edit("The genome-wide result is consistent with previous evidence", "sensitive to assumptions about sample overlap and the LD reference.", "sensitive to the assumed sample overlap and the multiple-testing family.")
edit("Supplementary Figure S1 Single-cell contextualisation", "Single-cell contextualisation of genes significant in both trait-specific MAGMA analyses", "Single-cell context of selected genes significant in both MAGMA analyses.")
edit("MAGMA gene-level analysis identified 147", "MAGMA genes—LDLRAD2, WNT4, GREB1, ESR1, SYNE1, VIP and ARL14EP—mapped", "MAGMA genes (LDLRAD2, WNT4, GREB1, ESR1, SYNE1, VIP and ARL14EP) mapped")
edit("Summary statistics were harmonised by rsID", "are reported above and in Supplementary Table S1.", "are reported in Table 1 and Supplementary Table S1.")
'''
a='# ---------------------------------------------------------------- figures 3 and 4'
r(a,extra+'\n'+a)
b='# ---------------------------------------------------------------- gene symbols in italics'
r(b,'''for h in ["Gene-level association", "MAGMA identifies genes significant in both disorders", "Single-cell expression and heritability enrichment point", "Fine-mapping compares signals between the disorders at selected regions", "Integrated genetic analyses identify regional signals"]:
    for p in doc.paragraphs:
        if p.text.strip().startswith(h) and len(p.text) < 130:
            for r_ in p.runs:
                r_.bold = True
for p in list(doc.paragraphs):
    if p.text and not p.text.strip() and "\\t" in p.text:
        p._p.getparent().remove(p._p)
'''+b)
c='hl = docx.Document(DST)'
r(c,'''r1 = docx.Document(FIG + "Revsions_round1_submitted/Manuscript_R1_clean.docx")
old_r1 = {norm(p.text) for p in r1.paragraphs}
FORCE = ("Vijayachitra Modhukur1", "OÜ Genorama", "Celvia", "Department of Obstetrics", "Nalam", "The present study was funded", "Supplementary Figure S1", "Supplementary Table S1")
'''+c)
d='''        for r in p.runs:
            r.font.highlight_color = WD_COLOR_INDEX.YELLOW
        n += 1'''
r(d,d+'''
    elif norm(p.text) and p.text.strip().startswith(FORCE) and norm(p.text) not in old_r1:
        for r in p.runs:
            r.font.highlight_color = WD_COLOR_INDEX.YELLOW
        n += 1''')
open('edit_manuscript_R2.py','w').write(s)
print('patched')
