"""Apply the second-round (R2) revisions to a COPY of the 30 Oct 2026 manuscript.

Source : Funciton_int_ge/Manuscript_R1_clean_affiliations_SekMo_2026-10-30.docx   (not modified)
Output : Funciton_int_ge/Manuscript_R2_revised.docx and Manuscript_R2_changes_highlighted.docx

Each edit is located by a unique text anchor; the script stops if an anchor is missing or ambiguous.
"""
import copy, re, shutil
import docx
from docx.enum.text import WD_COLOR_INDEX
from docx.oxml.ns import qn
from PIL import Image

FIG = "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/"
SRC = FIG + "Manuscript_R1_clean_affiliations_SekMo_2026-10-30.docx"
DST = FIG + "Manuscript_R2_revised.docx"
shutil.copy2(SRC, DST)
doc = docx.Document(DST)
log = []


def find(prefix):
    hits = [p for p in doc.paragraphs if p.text.strip().startswith(prefix)]
    assert len(hits) == 1, f"{len(hits)} paragraphs start with {prefix!r}"
    return hits[0]


def sub(p, old, new):
    text = "".join(r.text for r in p.runs)
    assert text.count(old) == 1, f"{text.count(old)} matches of {old[:60]!r} in {text[:40]!r}"
    s = text.index(old); e = s + len(old); pos = 0; done = False
    for r in p.runs:
        a, b = pos, pos + len(r.text); pos = b
        if b <= s or a >= e:
            continue
        lo, hi = max(s, a) - a, min(e, b) - a
        r.text = r.text[:lo] + ("" if done else new) + r.text[hi:]
        done = True


def edit(prefix, old, new):
    sub(find(prefix), old, new)
    log.append(prefix[:40])


# ---------------------------------------------------------------- Abstract
edit("Results: The genome-wide",
     "remained significant under all four overlap-aware sampling correlations and whenever estimable without overlap modelling.",
     "remained significant under all four overlap-aware sampling correlations (four regions under fixed correction families); three could not be tested without overlap.")
edit("Background: Uterine fibroids", "Whether this sharing is concentrated in particular genomic regions, genes and cell types is unclear.", "Whether it is concentrated in particular regions, genes and cell types is unclear.")
edit("Results: The genome-wide",
     "A stromal marker annotation was enriched for heritability in both disorders.",
     "A stromal marker annotation was enriched for heritability in both disorders, but the conditional result was marker-dependent.")
edit("Results: The genome-wide",
     "At ESR1/SYNE1, both disorders shared a single-variant credible set.",
     "At ESR1/SYNE1, both disorders shared a single-variant credible set (descriptive; no formal colocalisation).")

edit("Results: The genome-wide", "Of 88 testable genomic blocks", "Of 88 testable blocks")
edit("Conclusions: The two disorders show", "with the strongest detectable local evidence", "with the strongest local evidence")
edit("Conclusions: The two disorders show", "All genetic association and fine-mapping results derive from FinnGen Release 9 and require replication in independent cohorts.", "All genetic results derive from FinnGen Release 9 and require independent replication.")
edit("Uterine fibroids (uterine leiomyoma) and endometriosis", "Uterine fibroids (uterine leiomyoma)", "Uterine fibroids (uterine leiomyomas)")
edit("We conducted a comparative genetic epidemiologic study", "genetic epidemiologic study", "genetic epidemiological study")
edit("Methods: We analysed", "and Sum of Single Effects (SuSiE) fine-mapping was compared between the disorders.", "and Sum of Single Effects (SuSiE) fine-mapping compared the disorders.")
edit("Results: The genome-wide", "three could not be tested without overlap.", "three were untestable without overlap.")
edit("Results: The genome-wide", "Of 88 testable blocks", "Of 88 blocks with valid P values")
edit("Results: The genome-wide", "after Bonferroni correction, ten of them in both.", "after Bonferroni correction, ten in both.")
edit("Results: The genome-wide", "16 showed significant local genetic correlation after", "16 showed significant local correlation after")
edit("Results: The genome-wide", "three were untestable without overlap.", "three were untestable without overlap.".replace("three were untestable", "three untestable"))
# ---------------------------------------------------------------- Figure 1 caption
edit("Figure 1: Integrated", "Figure 1: Integrated", "Figure 1. Integrated")
edit("Figure 1. Integrated", "of which eight (five regions) were robust to the assumed sample overlap",
     "of which eight (five regions) remained significant under all overlap-aware sampling correlations when FDR was calculated within each run")

# ---------------------------------------------------------------- Methods
edit("FinnGen Release 9 GWAS summary statistics were obtained",
     "MAGMA used the munged-file sample sizes, and repeating it with the effective sample sizes gave the same significant genes (Supplementary Tables S1 and S2).",
     "MAGMA used the effective sample sizes (see Gene-level association).")
edit("Local genetic correlation was estimated with LAVA",
     "within the range of values tested in the sensitivity analysis.",
     "within the range of values tested in the sensitivity analysis. The cross-trait intercept reflects the covariance of the Z-score errors induced by participants shared between the two GWAS: shared cases and shared controls raise it, and participants who are cases in one study and controls in the other lower it. Dividing it by the square root of the product of the two intercepts gives the correlation of the sampling errors that LAVA requires. The intercept also absorbs other shared confounding, so the value is a model-based approximation of the overlap and not a count. Its estimation error is small (standard error of the standardised value about 0.008) relative to the differences between the values examined (Supplementary Table S5, Panel G).")
edit("Local genetic correlation was estimated with LAVA",
     "Blocks significant at FDR < 0.05 under every non-zero value, and whenever estimable without overlap, were classified as stable. Other significant blocks were classified as sensitivity-dependent, and adjacent significant blocks were interpreted as a single region.",
     "FDR was calculated within each run. Blocks significant at FDR < 0.05 under every overlap-aware value, and also without overlap where a valid test was returned, were classified as stable, so the classification refers to the tests that could be performed. Because the number of valid tests changed with the specification (26, 88, 93, 61 and 59), the comparison was repeated with two fixed correction families (Benjamini–Hochberg over the 88 primary tests, and over all 148 candidate rows with rows lacking a valid P value assigned P = 1, the most conservative option) and within the 14 blocks that were estimable under all five specifications. Other significant blocks were classified as sensitivity-dependent, and adjacent significant blocks were interpreted as a single region.")
edit("Local genetic correlation was estimated with LAVA",
     "Complete results are provided in Supplementary Tables S1 and S5.",
     "Complete results, including the raw P values, estimates and missing-result flags for all 148 rows in each specification, are provided in Supplementary Tables S1 and S5.")
edit("Gene-level association analysis was performed",
     "MAGMA was also run with the Release 9 effective sample sizes.",
     "The primary analysis, which supplies the gene statistics in Figure 3 and Supplementary Table S2, used the Release 9 effective sample sizes. An earlier run with the constant sample sizes of the munged files is reported as a sensitivity analysis (Supplementary Tables S1 and S2).")
edit("Publicly available single-cell RNA-seq data were used",
     "The coefficient test was one-sided because the hypothesis, as in that framework, was prespecified as positive enrichment; two-sided values are also reported (Supplementary Table S4).",
     "The coefficient test was one-sided, testing positive enrichment as in that framework; we have no record that this direction was fixed before the results were examined, and it was not specified in a registered analysis plan, so two-sided values are also reported (Supplementary Table S4).")
edit("To assess uncertainty in the GSE218044 labels",
     "and S-LDSC was repeated with these annotations.",
     "and S-LDSC was repeated with these annotations. The marker lists of both sets and the code for the cell-level reannotation are provided in Supplementary Table S4 (Panel J) and the code repository.")
edit("Fine-mapping used two sources",
     "only high-purity credible sets (FinnGen good_cs) were counted.",
     "only high-purity credible sets (FinnGen good_cs, which in the downloaded files corresponds to a minimum pairwise r² of at least 0.25 within the set) were counted. Credible sets were read from the files finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.SUSIE.cred.summary.tsv and finngen_R9_N14_ENDOMETRIOSIS.SUSIE.cred.summary.tsv and the matching .snp.filter.tsv files.")
edit("Fine-mapping used two sources",
     "Credible sets were considered shared between disorders when they contained the same variants.",
     "Credible sets were considered shared between disorders when they contained the same variants. Because both endpoints come from the same FinnGen release and may share participants, and because a region can contain several signals that the two analyses may resolve differently, matching credible sets are reported as a descriptive observation: they do not quantify the evidence for a shared causal variant, and no formal colocalisation was performed.")

# ---------------------------------------------------------------- Results: LAVA
edit("The number of significant blocks depended",
     "Without overlap modelling, only 26 blocks could be tested, and all eight primary blocks among them remained significant.",
     "Without overlap modelling, only 26 blocks could be tested; eight of the 16 primary blocks were among them, and all eight were significant.")
edit("The number of significant blocks depended",
     "Eight blocks were significant under every overlap-aware value and whenever estimable without overlap. They formed five regions",
     "Eight blocks were significant under every overlap-aware value. Five of them (three ESR1/SYNE1 blocks, the WNT4 block and the chromosome 11 block at 32.1 Mb) were also significant without overlap modelling. The other three (the fourth ESR1/SYNE1 block, the GREB1 block and the chromosome 4 block) returned no valid estimate in that analysis and were therefore not tested successfully under that assumption. With a fixed correction family (Benjamini–Hochberg over the 88 primary tests), seven of the eight remained significant at every overlap-aware value; the chromosome 11 block did not (FDR 0.075 at 0.369). With the most conservative family (all 148 candidate rows, rows lacking a valid P value assigned P = 1), six remained significant, because the fourth ESR1/SYNE1 block was also lost (FDR 0.051 at 0.369). The core is therefore eight blocks in five regions when FDR is calculated within each run, and six or seven blocks in four regions (ESR1/SYNE1, WNT4, GREB1 and chromosome 4) under the fixed families. Only 14 blocks were estimable under all five specifications, and six of them were significant under every specification when corrected within that set: the three ESR1/SYNE1 blocks, WNT4 and the two chromosome 11 blocks (Supplementary Table S5, Panels D–F). The eight blocks formed five regions")
edit("The number of significant blocks depended",
     "rather than a precise ranking of local correlation.",
     "rather than a precise ranking of local correlation, and the five-region core is conditional on the specifications examined and on the tests that could be performed.")
edit("Figure 2. Local genetic correlation",
     "(stable core, eight blocks in five regions).",
     "(stable core with FDR calculated within each run, eight blocks in five regions; three of the eight returned no valid estimate without overlap modelling).")
edit("MAGMA gene-level analysis identified 147",
     "With the Release 9 effective sample sizes, MAGMA identified the same significant genes in both analyses, with gene Z statistics changing by at most 0.82 for fibroids and 0.05 for endometriosis (Supplementary Table S2).",
     "In an earlier run with the constant sample sizes of the munged files, MAGMA identified the same significant genes in both analyses, with gene Z statistics differing by at most 0.82 for fibroids and 0.05 for endometriosis (Supplementary Table S2).")
edit("Figure 3. Gene-level association", "Results are in Supplementary Table S2.",
     "Gene statistics are from the analysis with the Release 9 effective sample sizes. Results are in Supplementary Table S2.")

# ---------------------------------------------------------------- Results: S-LDSC
edit("Cell-type heritability enrichment analysis using S-LDSC",
     "(Figure 4; Supplementary Table S4). With two-sided tests",
     "(Figure 4; Supplementary Table S4). Figure 4B shows the one-sided coefficient tests with the original marker set. With two-sided tests")
edit("Cell-type heritability enrichment analysis using S-LDSC",
     "The stromal enrichment is therefore robust to label uncertainty, but the evidence that the stromal annotation contributes heritability beyond the baseline annotations is modest.",
     "Stromal enrichment therefore persisted under label uncertainty, whereas the conditional coefficient (baseline v1.2) depended on both the marker assignment and the test direction: with the original markers it was significant one-sided but not two-sided in either disorder, and with the re-derived markers it was not significant one-sided in either disorder (Supplementary Table S4, Panels H(iii) and H(iv)). Enrichment (the proportion of heritability relative to the proportion of SNPs) and the conditional coefficient (the contribution beyond the baseline annotations) address different questions and should not be read as interchangeable. Because most cells in the smooth-muscle cluster were relabelled as stromal, the stromal and smooth-muscle signals cannot be cleanly separated in this dataset, and the evidence that the stromal annotation contributes heritability beyond the baseline annotations is modest.")
edit("Figure 4. Cell-type heritability enrichment",
     "(B) Conditional coefficient Z score under the baseline v1.2 model; asterisks mark one-sided coefficient FDR < 0.05.",
     "(B) Conditional coefficient Z score under the baseline v1.2 model for the original marker set; asterisks mark one-sided coefficient FDR < 0.05. The coefficient result depends on the marker assignment and the test direction (stromal coefficient FDR with the original markers, two-sided: 0.089 for fibroids and 0.054 for endometriosis; with the re-derived markers, one-sided: 0.11 and 0.062), and enrichment (A) and the coefficient (B) test different quantities; see Supplementary Table S4, Panels H(iii) and H(iv).")

# ---------------------------------------------------------------- Results: fine-mapping and integration
edit("Fine-mapping distinguishes shared and distinct signals", "Fine-mapping distinguishes shared and distinct signals at selected regions",
     "Fine-mapping compares signals between the disorders at selected regions")
edit("SuSiE fine-mapping with in-sample LD, published by FinnGen",
     "rs58415480 formed a single-variant credible set in both disorders, consistent with a shared signal.",
     "rs58415480 formed a single-variant credible set in both disorders. This is a descriptive observation: no formal colocalisation or independent replication was performed, both endpoints come from the same release with possible sample overlap, and matching credible sets do not quantify the evidence for a shared causal variant.")
edit("Across the genetic analyses, uterine fibroids and endometriosis showed",
     "in-sample fine-mapping identified a shared single-variant credible set at ESR1/SYNE1.",
     "in-sample fine-mapping identified a shared single-variant credible set at ESR1/SYNE1, a descriptive observation without formal colocalisation.")

# ---------------------------------------------------------------- Discussion and Conclusions
edit("Across the genetic analyses, uterine fibroids and endometriosis showed", "eight of which formed a robust core across five genomic regions", "eight of which, across five genomic regions, remained significant under all overlap-aware sampling correlations when FDR was calculated within each run (six or seven blocks in four regions under fixed correction families)")
edit("This study characterised the genetic relationship", "were robust to assumptions about sample overlap", "remained significant under all overlap-aware sampling correlations with FDR calculated within each run")
edit("The genome-wide result is consistent with previous evidence",
     "Only eight, in five regions, remained significant across all the sampling correlations tested, and the other eight depended on the assumed sample overlap. The evidence therefore supports a core of shared regions, notably ESR1/SYNE1, WNT4 and GREB1, rather than 16 independent loci.",
     "Eight, in five regions, remained significant across all the overlap-aware sampling correlations, and the other eight depended on the assumed sample overlap. Three of the eight, including the GREB1 block, could not be tested without overlap modelling, chromosome 11 at 32.1 Mb lost significance under fixed correction families, and the fourth ESR1/SYNE1 block did so under the most conservative one. The evidence therefore supports a core of shared regions, most securely ESR1/SYNE1 and WNT4, with GREB1 supported only under the overlap-aware specifications, rather than 16 independent loci.")
edit("The cellular analyses provide biological context",
     "Cell-level reannotation also showed that stromal and smooth-muscle labels were not cleanly separable in this dataset.",
     "Cell-level reannotation also showed that stromal and smooth-muscle labels were not cleanly separable in this dataset (98% of smooth-muscle cluster cells were relabelled).")
edit("Fine-mapping helped to separate shared signals from nearby ones.",
     "Fine-mapping helped to separate shared signals from nearby ones.", "Fine-mapping compared signals between the disorders at selected regions.")
edit("Fine-mapping compared signals between the disorders",
     "consistent with a shared signal, although formal colocalisation is needed to establish a common causal variant.",
     "which we regard as descriptive: no formal colocalisation or independent replication was performed, both endpoints were analysed in the same FinnGen release with possible sample overlap, and matching credible sets do not quantify the evidence for a shared causal variant.")
edit("Uterine fibroids and endometriosis share a moderate component",
     "The strongest local evidence occurs at ESR1/SYNE1, WNT4 and GREB1, alongside",
     "The strongest local evidence occurs at ESR1/SYNE1 and WNT4, with GREB1 supported under the overlap-aware specifications, alongside")
edit("Uterine fibroids and endometriosis share a moderate component",
     "in-sample fine-mapping identified a shared single-variant credible set at ESR1/SYNE1.",
     "in-sample fine-mapping identified a shared single-variant credible set at ESR1/SYNE1, a descriptive observation without formal colocalisation or independent replication.")

# ---------------------------------------------------------------- Availability, supplementary list
edit("Analysis code and project documentation",
     "is archived at Zenodo (https://doi.org/10.5281/zenodo.23692773). ",
     "is archived at Zenodo (https://doi.org/10.5281/zenodo.23912805). The second-round analyses (LAVA results for every specification, the comparison of marker sets and test directions in S-LDSC, and the MAGMA run with effective sample sizes) are included in the same release. ")
edit("Analysis code and project documentation", "the version used for this article (release v1.0-R1.5)", "the version used for this article (release v1.0-R2)")
edit("Supplementary Table S2 Complete MAGMA",
     "effective-sample-size sensitivity results and the LAVA blocks overlapping the ten shared genes.",
     "the first-run results with the munged-file sample sizes (sensitivity) and the LAVA blocks overlapping the ten shared genes. The primary gene statistics use the Release 9 effective sample sizes.")
edit("Supplementary Table S4 Single-cell context",
     "S-LDSC with an all-genes control annotation.",
     "S-LDSC with an all-genes control annotation, the dependence of the conditional coefficient on marker set and test direction, and the marker gene lists.")
edit("Supplementary Table S5 LAVA results",
     "and the blocks significant only without overlap modelling.",
     "the blocks significant only without overlap modelling, all 148 bivariate rows in each specification with raw P values and missing-result flags, the comparison under a fixed correction family and a common estimable set, and the derivation of the sampling correlations.")


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

# ---------------------------------------------------------------- figures 3 and 4
Image.MAX_IMAGE_PIXELS = None
NEW = {"2": FIG + "revision_figures/Figure_2_LAVA_R2.png", "3": FIG + "revision_figures/Figure_3_MAGMA_R2.png", "4": FIG + "revision_figures/Figure_4_SLDSC_R2.png"}
imgs = [p for p in doc.paragraphs if p._p.xpath(".//a:blip")]
assert len(imgs) == 4
for key, para in zip(("1", "2", "3", "4"), imgs):
    if key not in NEW:
        continue
    for b in para._p.iter(qn("a:blip")):
        part = doc.part.related_parts[b.get(qn("r:embed"))]
        assert str(part.partname).endswith(".png"), part.partname
        part._blob = open(NEW[key], "rb").read()
        w, h = Image.open(NEW[key]).size
        for ext in para._p.iter(qn("wp:extent"), qn("a:ext")):
            ext.set("cx", "5727700"); ext.set("cy", str(int(5727700 * h / w)))

for h in ["Gene-level association", "MAGMA identifies genes significant in both disorders", "Single-cell expression and heritability enrichment point", "Fine-mapping compares signals between the disorders at selected regions", "Integrated genetic analyses identify regional signals"]:
    for p in doc.paragraphs:
        if p.text.strip().startswith(h) and len(p.text) < 130:
            for r_ in p.runs:
                r_.bold = True
for p in list(doc.paragraphs):
    if p.text and not p.text.strip() and "\t" in p.text:
        p._p.getparent().remove(p._p)
# keep each figure caption on one page with its figure, and widen the first column of Table 1
for _i, _p in enumerate(doc.paragraphs):
    if _p._p.xpath(".//pic:pic"):
        _p.paragraph_format.keep_with_next = True
        doc.paragraphs[_i + 1].paragraph_format.keep_together = True
_t = doc.tables[0]
_w = [2800, 3153, 1474, 1587]
for _g, _x in zip(_t._tbl.tblGrid.findall(qn("w:gridCol")), _w):
    _g.set(qn("w:w"), str(_x))
for _row in _t.rows:
    for _c, _x in zip(_row.cells, _w):
        _c.width = docx.shared.Emu(_x * 635)
# ---------------------------------------------------------------- gene symbols in italics
GENES = ["ESR1", "SYNE1", "WNT4", "GREB1", "WT1", "DNM3", "EFR3B", "VIP", "LDLRAD2", "PLEKHM3", "ARL14EP"]
GENE_RE = re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(GENES) + r")(?![A-Za-z0-9_])")


def italicise(para):
    for r in list(para.runs):
        if r.italic or not GENE_RE.search(r.text):
            continue
        parts = GENE_RE.split(r.text)
        r.text = parts[0]
        prev = r
        for k, piece in enumerate(parts[1:], 1):
            if not piece:
                continue
            nr = copy.deepcopy(r._r)
            prev._r.addnext(nr)
            run = docx.text.run.Run(nr, para)
            run.text = piece
            run.italic = True if k % 2 == 1 else r.italic
            prev = run


end = find("References")
for para in doc.paragraphs:
    if para._p is end._p:
        break
    italicise(para)
doc.save(DST)

# ---------------------------------------------------------------- highlighted copy (baseline = source)
norm = lambda t: re.sub(r"\s+", " ", t).strip()
base = docx.Document(SRC)
old = {norm(p.text) for p in base.paragraphs}
r1 = docx.Document(FIG + "Revsions_round1_submitted/Manuscript_R1_clean.docx")
old_r1 = {norm(p.text) for p in r1.paragraphs}
FORCE = ("Vijayachitra Modhukur1", "OÜ Genorama", "Celvia", "Department of Obstetrics", "Nalam", "The present study was funded", "Supplementary Figure S1", "Supplementary Table S1")
hl = docx.Document(DST)
n = 0
for p in hl.paragraphs:
    if norm(p.text) and norm(p.text) not in old:
        for r in p.runs:
            r.font.highlight_color = WD_COLOR_INDEX.YELLOW
        n += 1
    elif norm(p.text) and p.text.strip().startswith(FORCE) and (norm(p.text) not in old_r1 or p.text.startswith('Nalam')):
        for r in p.runs:
            r.font.highlight_color = WD_COLOR_INDEX.YELLOW
        n += 1
hl.save(FIG + "Manuscript_R2_changes_highlighted.docx")
print(len(log), "edits; highlighted paragraphs:", n)
