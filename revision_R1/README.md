# Revision R1: analyses for the response to reviewers

Code and derived results added for the first revision of:

> Modhukur V, Lingasamy P, Patel N, Salumets A. *Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis.* Functional & Integrative Genomics (under revision, 2026).

All analyses use public summary-level data. No individual-level or controlled-access data are included.

## Inputs

| Input | Identifier | Source |
|---|---|---|
| Uterine fibroid GWAS | FinnGen R9 `CD2_BENIGN_LEIOMYOMA_UTERI` (31,661 cases, 179,209 controls), file `finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz` | https://finngen.gitbook.io/documentation/r9/data-download |
| Endometriosis GWAS | FinnGen R9 `N14_ENDOMETRIOSIS` (15,088 cases, 107,564 controls), file `finngen_R9_N14_ENDOMETRIOSIS.gz` | as above |
| FinnGen R9 SuSiE fine-mapping (in-sample LD) | `finngen_R9_<endpoint>.SUSIE.cred.summary.tsv`, `.SUSIE.snp.filter.tsv` | `scripts/finemap_finngen/download_finngen_r9_susie.sh` |
| LD reference | 1000 Genomes Phase 3 European PLINK files | 1000 Genomes Project |
| S-LDSC baselines | baseline v1.2, baseline-LD v2.2 | Zenodo record 7768714 |
| Single-cell markers | GSE218044 (uterine tissue, four individuals with adenomyosis), processed in scEndoExplorer | NCBI GEO; Chen et al., *Front Immunol* 2022;13:945504 |
| Fibroid single-cell context | GSE162122 (fibroid and matched myometrium) | NCBI GEO; Goad et al., *Hum Reprod* 2022 |

Coordinates: FinnGen statistics and fine-mapping are GRCh38; the 1000 Genomes reference, LAVA blocks and MAGMA gene boundaries are GRCh37. Variants were matched by rsID and alleles.

Scripts contain the absolute paths of the analysis machine; adjust them to your local copies.

## Output-to-script map

| Output (results/) | Script (scripts/) | Supplementary table / figure |
|---|---|---|
| `ldsc/rg_{submitted,total,neff}.log` | `ldsc/ldsc_rg_sample_size_sensitivity.sh` | Table S1 Panel B |
| `lava/LAVA_corrected_N_ov01724_all_*.csv`, `lava/res_p1_c*.csv`, `lava/qc_p1_c*.csv` | `lava/run_lava.R` (2,495 blocks, 4 chunks; `chain_c3.sh`) with `info_corrected.tsv`, `ov_0.1724.tsv` | Table S1 Panels C–D; Figure 2 |
| `lava/res_grid_corr_{0.0,0.2590,0.35}.csv`, `lava/LAVA_16_blocks_*.csv` | `lava/run_lava_subset.R` with `ov_*.tsv` | Table S5; Figure 2B |
| `lava/res_repro_legacy.csv` | `lava/run_lava_subset.R` with `info_legacy.tsv`, `ov_original_0.17240.tsv` | Table S5 Panel A (reproduction row) |
| `magma/{fibroids,endo}_neff.genes.out` | `magma/magma_effective_n.sh` | Table S2 |
| `magma/MAGMA_cross_trait_comparison.csv` | `supplementary/build_supplementary_R1.py` (Panel A) | Table S2; Figure 3 |
| `susie/SuSiE_L_sensitivity_rerun.csv`, `SuSiE_LD_concordance_diagnostics.csv`, `SuSiE_credible_set_membership_L10.csv`, `kriging_flagged_*.csv` | `susie/susie_rerun.R` | Table S3 Panels D–F; Supplementary Figure S2 |
| `finemap_finngen/FinnGen_insample_*.csv` | `finemap_finngen/compare_finngen_insample.py` | Table S3 Panels A–C |
| `sldsc/SLDSC_full_statistics_from_raw_results.csv`, `SLDSC_marker_gene_set_overlap.csv` | `supplementary/build_supplementary_R1.py` | Table S4 Panels B–D; Figure 4 |
| `reannotation/*` | `reannotation/reannotate_GSE218044.R` | Table S4 Panels E–H |
| `sldsc/SLDSC_original_vs_confident_markers.csv` | `sldsc/make_annot.py`, `sldsc/run_sldsc.sh`, `sldsc/chain_confident.sh` (see `sldsc/LDSC_PATCH.md`) | Table S4 Panel H(iii) |
| `build_map/*` | rsID lookup in the FinnGen R9 files and match to the 1000 Genomes `.bim` (inline commands); assembled into Table S6 by `supplementary/build_supplementary_R1.py` | Table S6 |
| Figures 2–4, Supplementary Figures S1–S2 | `figures/*` | Figures 2–4, S1, S2 |

## Key revised results

- **Genome-wide correlation:** rg = 0.511 (SE 0.069; GenomicSEM). Standard LDSC gives rg = 0.510 under submitted, total and effective sample sizes.
- **Local correlation (LAVA), with R9 case and control counts:**
  - 148 of 2,495 blocks reached the bivariate test, 88 of them with valid P values; 16 were FDR-significant.
  - These are the same 16 blocks as in the submitted analysis.
  - Nine blocks in five regions were significant under every overlap-aware sampling correlation (0.1724, 0.259, 0.35) and whenever estimable without overlap: ESR1/SYNE1, WNT4, GREB1, chr4 at 70.5 Mb, and chr11.
- **MAGMA:** 147 and 38 Bonferroni-significant genes; 10 significant in both analyses. Results are unchanged with effective sample sizes.
- **Fine-mapping (FinnGen in-sample SuSiE):**
  - The same credible variant appears in both disorders at ESR1/SYNE1 (rs58415480) and at chr11 (rs11031005/rs11031006).
  - GREB1 and WNT4 have non-overlapping credible sets.
  - SuSiE with the 1000 Genomes LD reference was inconsistent with the summary statistics (estimate_s_rss 0.56–0.90) and is a sensitivity analysis only.
- **S-LDSC:** the stromal marker annotation is enriched in both disorders.
  - With markers re-derived from cells whose per-cell and cluster labels agree, enrichment remains significant (FDR 0.012 and 0.040).
  - The conditional coefficient test is attenuated (FDR 0.11 and 0.062).
