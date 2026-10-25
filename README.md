# Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis

This repository contains analysis code, endpoint information, software details and permitted aggregate outputs supporting the manuscript:

> Modhukur V, Lingasamy P, Patel N, Salumets A. *Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis.* Functional & Integrative Genomics (under revision, 2026).

## Overview

Uterine fibroids and endometriosis are common, oestrogen-responsive gynaecological disorders with partly overlapping biological features. This study uses publicly available FinnGen Release 9 genome-wide association study (GWAS) summary statistics and public reference resources to characterise their shared and distinct common-variant genetic architecture.

The study is based entirely on summary-level and publicly available data. It does not use UK Biobank data, individual-level participant data or controlled-access biobank data.

The analytical framework includes:

1. Genome-wide genetic correlation using LD Score Regression (LDSC)
2. Local genetic correlation using LAVA, with sensitivity analyses for the assumed sample overlap
3. MAGMA gene-level association analysis
4. Fine-mapping: FinnGen's published Release 9 SuSiE results (in-sample LD) as primary, with SuSiE using the 1000 Genomes reference as a sensitivity analysis
5. Single-cell contextualisation using published datasets and stratified LDSC (S-LDSC) with single-cell marker annotations

The exploratory pathway annotation in the original submission was removed in revision.

## Data sources

### FinnGen GWAS summary statistics

- Uterine fibroids: `CD2_BENIGN_LEIOMYOMA_UTERI`
- Endometriosis: `N14_ENDOMETRIOSIS`
- Release: FinnGen R9
- Access: <https://finngen.gitbook.io/documentation/r9/data-download>

Endpoint metadata used in the manuscript report:

| Endpoint | Cases | Controls |
|---|---:|---:|
| Uterine fibroids | 31,661 | 179,209 |
| Endometriosis | 15,088 | 107,564 |

Analysis-specific sample-size parameters retained in archived scripts are documented separately from endpoint catalogue denominators and should not be treated as interchangeable.

### Linkage-disequilibrium reference

The 1000 Genomes Project Phase 3 European reference panel was used for linkage-disequilibrium-based analyses where specified in the manuscript.

### Single-cell resources

Published single-cell RNA-seq resources were used for expression contextualisation and annotation construction, including:

- GSE162122: uterine fibroid and myometrial tissue
- GSE218044: uterine tissue from individuals with adenomyosis (source of the S-LDSC marker annotations; processed in scEndoExplorer)

Published cell labels were used for the primary analysis. In revision, GSE218044 cells were also reannotated per cell (SingleR) as a sensitivity analysis; see `revision_R1/`.

## Data availability and access restrictions

This repository is intended to contain only:

- analysis code;
- software and environment information;
- public endpoint identifiers and data-source information;
- non-disclosive aggregate results;
- manuscript tables and figures where redistribution is permitted.

No individual-level participant records, participant identifiers, controlled-access biobank files or restricted genotype data are included or required for this Paper 1 workflow.

External datasets remain subject to the terms imposed by their original providers. Users should obtain all required source files through the original repositories and comply with the applicable data-use and citation requirements.

## Repository structure

```text
revision_R1/                   Code and derived results for the revised manuscript (start here)
  README.md                    Inputs, output-to-script map and key revised results
  scripts/                     LDSC, LAVA, MAGMA, fine-mapping, S-LDSC, reannotation, figures, supplementary tables
  results/                     Derived results underlying Supplementary Tables S1–S6 and the figures

scripts/                       Scripts from the original submission (LDSC, LAVA, MAGMA, single-cell context, figures)
figures/                       Figures as originally submitted
environment/                   R and Python package information
```

The revised analyses, figures and Supplementary Tables S1–S6 are produced by the scripts in [`revision_R1/`](revision_R1/README.md), which also lists the script behind each output.

## Main manuscript outputs (revised)

### Main figures

1. Figure 1: Integrated analytical framework
2. Figure 2: LAVA local genetic correlation and sensitivity to sample overlap
3. Figure 3: Gene-level architecture of uterine fibroids and endometriosis
4. Figure 4: Cell-type heritability enrichment (S-LDSC)

### Main table

1. Table 1: Public data sources used in the study

### Supplementary figures

1. Supplementary Figure S1: Single-cell context of genes significant in both MAGMA analyses
2. Supplementary Figure S2: SuSiE fine-mapping with the 1000 Genomes LD reference (sensitivity analysis)

### Supplementary tables

1. Supplementary Table S1: Sample sizes, LDSC and LAVA
2. Supplementary Table S2: MAGMA gene-level results
3. Supplementary Table S3: Fine-mapping (FinnGen in-sample SuSiE and 1000 Genomes-LD sensitivity analysis)
4. Supplementary Table S4: S-LDSC and single-cell context
5. Supplementary Table S5: LAVA sensitivity to the sampling correlation
6. Supplementary Table S6: Genome builds

## Software

| Tool | Version or implementation | Primary reference |
|---|---|---|
| LDSC | LDSC/GenomicSEM workflow | Bulik-Sullivan et al., *Nature Genetics* (2015) |
| LAVA | v0.1.5 | Werme et al., *Nature Genetics* (2022) |
| MAGMA | v1.10 | de Leeuw et al., *PLoS Computational Biology* (2015) |
| susieR | Sum of Single Effects fine-mapping | Wang et al., *Journal of the Royal Statistical Society: Series B* (2020) |
| S-LDSC | Stratified LD Score Regression | Finucane et al., *Nature Genetics* (2015) |
| R | See `environment/R_packages.txt` | — |
| Python | See `environment/requirements.txt`, if used | — |

Package versions, parameters, reference files and analysis-specific sample-size values should be preserved in the final publication archive.

## Reproducibility notes

- The analyses use public GWAS summary statistics and public reference resources; no controlled-access data are required.
- Sample sizes, harmonisation and sampling-correlation specifications are listed in Supplementary Table S1 and `revision_R1/README.md`.
- LAVA results are reproduced with the primary LDSC-standardised sampling correlation (0.1724); the values 0, 0.259, 0.35 and 0.369 are sensitivity analyses.
- Primary fine-mapping uses FinnGen's published Release 9 SuSiE credible sets (high-purity sets within ±500 kb of the regional lead variant). SuSiE with 1000 Genomes European LD was inconsistent with the summary statistics and is a sensitivity analysis only.
- Fine-mapping and cellular analyses are prioritisation or contextualisation layers and do not establish causal variants, genes or cell types.

## Citation

Until publication:

> Modhukur V, Lingasamy P, Patel N, Salumets A. *Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis.* Functional & Integrative Genomics (under revision, 2026).


## Correspondence

**Vijayachitra Modhukur**  
Department of Obstetrics and Gynaecology  
Institute of Clinical Medicine, University of Tartu, Estonia  
Email: vijayachitra.modhukur@ut.ee
