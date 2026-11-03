"""Point-by-point response to the second-round reviewer comments (R2). Yellow = placeholders to resolve before PDF export."""
import re
import docx
from docx.enum.text import WD_COLOR_INDEX
from docx.shared import Pt, Cm

OUT = "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/Response_to_Reviewers_FIG_R2.docx"
d = docx.Document()
d.styles["Normal"].font.name = "Calibri"
d.styles["Normal"].font.size = Pt(11)
for s in d.sections:
    s.left_margin = s.right_margin = Cm(2.3)


def runs(p, text):
    for part in re.split(r"(\[\[.*?\]\])", text):
        if not part:
            continue
        if part.startswith("[["):
            p.add_run(part[2:-2]).font.highlight_color = WD_COLOR_INDEX.YELLOW
        else:
            p.add_run(part)


def para(text, bold_prefix=None):
    p = d.add_paragraph()
    if bold_prefix:
        p.add_run(bold_prefix).bold = True
    runs(p, text)
    return p


def comment(text):
    p = d.add_paragraph()
    p.paragraph_format.left_indent = Cm(0.5)
    r = p.add_run(text)
    r.italic = True


def response(*paras):
    for i, t in enumerate(paras):
        para(t, "Response: " if i == 0 else None)


d.add_heading("Response to the Editor and Reviewers (second revision)", 0)
para("Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis", "Manuscript: ")
para("Functional & Integrative Genomics    Submission ID: 4b617de9-5b8f-4cdb-b503-85d3672d4058", "Journal: ")
para("Reviewer comments are reproduced in italics, followed by our response. Locations refer to the clean revised manuscript and Supplementary Tables S1 to S6.")

d.add_heading("Response to the Editor", 1)
para("Dear Editor,")
para("Thank you and both reviewers for the further careful review. We have revised the manuscript, the supplementary tables and Figures 3 and 4. The main changes are listed below.")
for t in [
    "LAVA: the supplement now gives the raw P value, estimate, confidence interval, FDR and reason for a missing P value for all 148 rows in each of the five specifications (Supplementary Table S5, Panel D). The comparison was repeated with a fixed correction family and within the 14 blocks estimable in all five specifications (Panel F), and the derivation and uncertainty of the sampling correlation are given in Panel G. Three of the eight stable blocks (the fourth ESR1/SYNE1 block, GREB1 and chromosome 4) have no estimate without overlap modelling and are now described as not tested successfully under that assumption. The comparison was repeated with two fixed correction families. Under the 88-test family the chromosome 11 block (32.1 Mb) does not remain significant, and under the most conservative family (all 148 candidate rows, rows without a valid P assigned P = 1) the fourth ESR1/SYNE1 block is also lost, so the core is eight blocks in five regions with FDR within each run and six or seven blocks in four regions with fixed families. The core is stated to be conditional on the specifications and tests available.",
    "Stromal S-LDSC: the text, Figure 4 and its legend now state that the conditional coefficient depends on the marker assignment and the test direction, and distinguish enrichment from the conditional coefficient. We removed the statement that the one-sided direction was prespecified. Marker lists and reannotation code are provided.",
    "Fine-mapping: the shared single-variant credible set at ESR1/SYNE1 is described as a descriptive observation, without formal colocalisation or independent replication, in the Abstract, Methods, Results, Discussion and Conclusions. The source files and selection rule are stated.",
    "MAGMA: the gene statistics in Figure 3 and Supplementary Table S2 now come from the run with the Release 9 effective sample sizes, consistently in the Methods, Results, figure legend and supplement. The earlier run is retained as a sensitivity analysis.",
    "Reviewer 2: the main text, Table 1, figure legends and supplementary files were checked for typographical, formatting and grammatical errors, and corrected (listed under Reviewer 2).",
]:
    d.add_paragraph(t, style="List Bullet")
para("Yours sincerely,\nVijayachitra Modhukur, on behalf of all authors")

d.add_heading("Reviewer 1", 1)
comment("The authors have addressed many of my concerns. They corrected the sample-size inputs, replaced the unreliable out-of-sample fine-mapping as the primary analysis, removed the pathway panels and clarified the tissue source of the cellular annotations. The revised supplement also makes the overlap sensitivity and annotation analyses easier to examine. A few points still affect the interpretation of the local and cell-type findings.")
response("We thank the reviewer for these further comments, which sharpened the LAVA, S-LDSC and fine-mapping interpretation.")

d.add_heading("Comment 1", 2)
comment("1. The LAVA sensitivity analysis is more informative, but the evidence for a stable set depends on both the assumed sampling correlation and which tests produce valid P values. Supplementary Table S5 reports 26, 88, 93, 61 and 59 valid tests across the five specifications, and each FDR is calculated within that run. Three of the eight blocks classified as stable have no estimate in the no-overlap analysis. Please provide the raw P values, estimates and failure reasons for all 148 bivariate rows under each specification, and repeat the comparison using a prespecified common estimable set and a fixed correction family where feasible. Make clear that the three blocks without a no-overlap estimate were not tested successfully under that assumption. The primary sampling correlation is obtained by standardising the LDSC intercept matrix, while the count-derived values depend on an unverified overlap count. Please explain why the standardised intercept is the appropriate sampling-error correlation for these binary phenotypes, including how overlapping cases and controls enter, and discuss uncertainty in that quantity. The five-region core should remain conditional on the examined specifications and the available tests.")
response(
    "We agree. All of the requested material is now in Supplementary Table S5.",
    "Raw results for all 148 rows. Panel D lists, for each of the five specifications (740 rows), the local heritabilities, the univariate P values, rho with its 95% interval, the raw P value, the FDR within the run and the FDR in a fixed family. Where no P value was returned, the row says so. LAVA returns the bivariate estimate as missing when it falls outside the parameter bound (|rho| > 1.25) or cannot be estimated and does not report which of the two occurred, so a more specific failure reason cannot be given. The two candidate blocks that failed during locus processing in every run (the covariance matrix was not positive definite) are listed in Panel E.",
    "Fixed correction family and common estimable set. Only 14 blocks are estimable in all five specifications, and eight of the 16 primary blocks lie outside that set, so the common set cannot support the whole core. We therefore report four views (Panel F): FDR within each run (as before); a fixed family of the 88 primary tests (Benjamini–Hochberg, applied to each specification); the most conservative fixed family, in which all 148 candidate rows enter and rows without a valid P value are assigned P = 1; and the 14-block common set. Under the 88-test family, seven of the eight stable blocks remain significant at every overlap-aware value; the chromosome 11 block at 32.1 Mb does not (FDR 0.075 at 0.369). Under the 148-row family six remain, because the fourth ESR1/SYNE1 block is also lost (FDR 0.051 at 0.369). Both fixed families leave four regions (ESR1/SYNE1, WNT4, GREB1 and chromosome 4). The Abstract, Methods, Results and Discussion state which definition supports the eight-block, five-region result (FDR within each run) and report the fixed-family results. Within the common set, the three ESR1/SYNE1 blocks, WNT4 and both chromosome 11 blocks are significant in every specification. These results are in Methods (Local genetic correlation) and Results.",
    "Blocks without a no-overlap estimate. The fourth ESR1/SYNE1 block (block_1077), the GREB1 block (block_0216) and the chromosome 4 block (block_0666) returned no valid estimate in the no-overlap analysis. The Results now say they were not tested successfully under that assumption, so their stability refers to the overlap-aware values only, and the Discussion and Conclusions now say that GREB1 is supported only under the overlap-aware specifications; ESR1/SYNE1 and WNT4 are the most securely supported regions. The Abstract states that three of the eight could not be tested without overlap.",
    "Why the standardised intercept. LAVA requires the correlation between the sampling errors of the two sets of SNP effect estimates. The cross-trait LDSC intercept estimates the covariance of the Z-score errors caused by shared participants; shared cases and shared controls raise it, and participants who are cases in one study and controls in the other lower it. Dividing by the square root of the product of the two single-trait intercepts converts that covariance into a correlation, 0.1724. Both phenotypes are binary and case-control based, so the intercept reflects the net effect of case-case, control-control and case-control overlap and does not require the number of shared participants. The intercept also absorbs any shared confounding, so the value is a model-based approximation of the overlap, not a count.",
    "Uncertainty. A delta-method approximation from the intercept standard errors gives a standard error of about 0.008 (95% interval about 0.157 to 0.188), ignoring covariances between the intercept estimates; the fit with effective sample sizes gives 0.1716. This estimation error is small relative to the differences between the values examined (0 to 0.369). The larger uncertainty is whether the intercept represents the overlap, which the count-derived values address only partly because the overlap count (59,264) could not be verified. This is now stated in Methods and in Supplementary Table S5, Panel G.",
    "Conditional core. The Results now state that the five-region core is conditional on the specifications examined and on the tests that could be performed.",
)

d.add_heading("Comment 2", 2)
comment("2. The reannotation analysis materially changes the strength of the conditional stromal result. In Supplementary Table S4, panel H, the re-derived stromal annotation remains enriched, whereas its one-sided conditional coefficient FDR is 0.11 for fibroids and 0.062 for endometriosis under baseline v1.2. Under the original annotations, the fibroid two-sided coefficient FDR is 0.089, compared with 0.045 for the one-sided test. Figure 4B displays only the original one-sided results. Please make the dependence on marker assignment and test direction clear in the figure legend and Results, alongside the distinction between enrichment and the conditional coefficient. The 98% disagreement for the smooth-muscle cluster also warrants caution when distinguishing stromal from smooth-muscle signals. Explain whether the one-sided direction was fixed before examining these results and provide the marker lists and cell-level analysis code needed to reproduce the sensitivity comparison.")
response(
    "We agree with the reviewer's reading of Supplementary Table S4, Panel H.",
    "Figure and Results. The legend of Figure 4 now states that panel B shows the one-sided coefficient tests with the original marker set, gives the two-sided FDR for the original markers (0.090 and 0.054) and the one-sided FDR for the re-derived markers (0.11 and 0.062), and states that enrichment (panel A) and the conditional coefficient (panel B) test different quantities. The panel titles and footnote of the figure were changed accordingly. The Results now state that stromal enrichment persisted with the re-derived markers, whereas the conditional coefficient depended on the marker assignment and the test direction, and that the evidence for a contribution beyond the baseline annotations is modest. A new Panel H(iv) of Supplementary Table S4 gives enrichment and the one- and two-sided coefficient P values and FDRs for both marker sets and both baseline models. The Results also explain that two-sided FDRs can be smaller than one-sided ones because the two-sided family counts annotations with strongly negative coefficients as rejections; with the re-derived markers the endometriosis stromal coefficient FDR is 0.041 two-sided and 0.062 one-sided (fibroids 0.21 and 0.11). We report this and do not interpret it as stronger evidence.",
    "Smooth muscle. The Results and Discussion state that, because most cells in the smooth-muscle cluster were relabelled as stromal (98% received a different label), the two signals cannot be cleanly separated in this dataset.",
    "Direction of the test. We cannot document that the one-sided direction was fixed before the results were examined: the first one-sided coefficient P values in our analysis files date from the October revision, after the coefficient z-scores were available, and no analysis plan was registered. We used the one-sided test because it is the convention of the framework of Finucane et al. for cell-type enrichment. Our earlier response described it as prespecified, which we cannot support. We have removed that word from the Methods, the figure legend and Supplementary Table S4, now state that there is no record that the direction was fixed beforehand, and report two-sided values alongside (Supplementary Table S4, Panel H(iv)).",
    "Reproducibility. The marker lists of both sets (top 500 per cell class, with rank, log2 fold change and adjusted P) are in Supplementary Table S4, Panel J. The cell-level reannotation code, the per-cell labels and the S-LDSC scripts are in the archived code release (see the Availability statement).",
)

d.add_heading("Comment 3", 2)
comment("3. FinnGen in-sample fine-mapping is a useful improvement, and the prior claim of allelic complexity has been withdrawn. The shared single-variant credible set at ESR1/SYNE1 is still highlighted in the Abstract and Conclusions, although both traits were analysed in FinnGen with possible sample overlap and matching credible sets alone do not quantify evidence for a shared causal variant. Please make the absence of formal colocalisation and independent replication clear next to this statement. Specify the source files, the precise credible-set selection rule and purity threshold, and explain how possible multiple signals and overlapping samples affect any cross-trait comparison. If colocalisation is not feasible, retain the shared-credible-set observation as descriptive.")
response(
    "We agree and have retained the observation as descriptive. Formal colocalisation is not feasible with the public data used here: we do not have the in-sample LD matrices behind the FinnGen credible sets, and the 1000 Genomes reference was inconsistent with the summary statistics in these regions.",
    "Wording. The Abstract now reads “both disorders shared a single-variant credible set (descriptive; no formal colocalisation)”, and the Conclusions state that this is a descriptive observation without formal colocalisation or independent replication. The Results and Discussion say the same and add that both endpoints come from the same release with possible sample overlap.",
    "Source files and rule. The Methods name the files finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.SUSIE.cred.summary.tsv and finngen_R9_N14_ENDOMETRIOSIS.SUSIE.cred.summary.tsv and the matching .snp.filter.tsv files. A credible set was counted when its lead variant lay within ±500 kb of the regional GWAS lead variant and the FinnGen field good_cs was true. In the downloaded files good_cs is true exactly when the minimum pairwise r² within the set is at least 0.25 (|r| ≥ 0.5); we inferred this from the 141 sets in the two files (smallest r² among good sets 0.251, largest among the others 0.243), and the Methods report it as the observed correspondence. The rule and files are also stated in Supplementary Table S3, Panel A.",
    "Multiple signals and overlap. The Methods now explain that a region can contain several signals that the two analyses resolve differently, that shared participants make the two sets of summary statistics, and so their credible sets, more alike than independent samples would, and that matching credible sets therefore do not quantify the evidence for a shared causal variant.",
)

d.add_heading("Minor comments", 2)
comment("1. The Methods state that MAGMA used the sample sizes from the earlier configuration, with Release 9 effective sample sizes as a sensitivity analysis. The response says MAGMA was rerun with the corrected values. Please state consistently which analysis supplies the primary gene statistics in Figure 3 and Supplementary Table S2, and identify the corresponding input and output files.")
response(
    "Thank you; the two statements were inconsistent. The primary gene statistics in Figure 3 and Supplementary Table S2 (Panel A, first columns) now come from the MAGMA run with the Release 9 effective sample sizes (107,629 for fibroids and 52,928 for endometriosis). The earlier run with the constant sample sizes of the munged files is retained as a sensitivity analysis (last four columns of Panel A). The Methods (two places), Results, Figure 3 legend, Supplementary Table S1 and the Supplementary Information list were changed to say this. The significant genes are identical (147, 38 and 10 in both), so no result changed. The MAGMA input, output files and scripts for the effective-sample-size run are in the archived code release.",
)
comment("2. Please identify the archived files containing the complete LAVA results and the exact release or commit used to produce each revised table and figure.")
response(
    "The complete LAVA results are in Supplementary Table S5, Panel D (all 148 rows in each of the five specifications) and as CSV files in the archived code release, whose README maps each table and figure to its input files, script and output files. The revised tables and figures were produced from release v1.0-R2 (Zenodo DOI 10.5281/zenodo.23912805); the Availability statement of the manuscript names the release.",
)

d.add_heading("Reviewer 2", 1)
comment("In this revised manuscript, Revision 2 (i.e., R2), the authors appear to have addressed all major comments and minor comments, and the authors should perform thorough checking to correct all possible errors in the main text of manuscript (including Table 1 and Figures 1-4 and all Supplementary Tables and all Supplementary Figures), including any possible typographical errors, formatting errors (for tables and figures), and grammatical errors.")
response(
    "Thank you. We re-read the main text, Table 1, the legends of Figures 1 to 4, Supplementary Figures S1 and S2 and the text of Supplementary Tables S1 to S6, and ran automated checks for repeated words, double spaces, spacing before punctuation, missing spaces and unrecognised words. We corrected the following:",
)
for t in [
    "The caption of Figure 1 used a colon (“Figure 1:”) whereas the other captions use a full stop; it now reads “Figure 1.”.",
    "The author and affiliation block of the Supplementary Figures file differed from the manuscript (affiliation numbering, and the Genorama affiliation was missing); it now matches the manuscript.",
    "The heading and first sentence of the fine-mapping sections no longer say that fine-mapping “distinguishes” or “separates” shared and distinct signals, which overstated the evidence; they now say that it compares signals.",
    "“Robust” was replaced by a precise statement (significant under all overlap-aware sampling correlations) where it described the LAVA core.",
    "Supplementary Table S1 and the Supplementary Information list now describe the MAGMA primary analysis consistently with the main text; the abstract remains within 250 words.",
    "The Introduction and Methods wording was corrected (“uterine leiomyomas”, “epidemiological”), and the Supplementary Figure S1 legend now describes colour intensity as mean log1p-transformed CP10K expression.",
    "In Supplementary Table S4 the panels are in sequence (the marker lists are Panel J), and the note to Panel J explains that adjusted P values shown as 0 are numerical underflow, not exact zeros.",
    "The caption of Figure 1 now states that the eight stable blocks refer to FDR calculated within each run.",
    "New table panels have complete descriptive titles, and all P values below 0.001 in the workbook are formatted in scientific notation.",
]:
    d.add_paragraph(t, style="List Bullet")
para("The legends of Figures 1 to 4 and Supplementary Figures S1 and S2 were checked; Figures 3 and 4 and Supplementary Figure S1 were redrawn (Figure 3 from the effective-sample-size run, with the overlapping gene labels in panel C separated; Figure 4 with clearer panel B labelling; Supplementary Figure S1 with the clipped labels in panel B no longer cut off).")
d.save(OUT)
print("wrote", OUT)
