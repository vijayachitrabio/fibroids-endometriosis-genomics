# Integrative genomic analyses reveal regional and cellular convergence between uterine fibroids and endometriosis

This repository contains analysis code, endpoint information, software details and permitted aggregate outputs supporting the manuscript:

> Modhukur V, Lingasamy P, Patel N, Salumets A. *Integrative genomic analyses reveal regional and cellular convergence between uterine fibroids and endometriosis.* Manuscript in preparation (2026).

## Overview

Uterine fibroids and endometriosis are common, oestrogen-responsive gynaecological disorders with partly overlapping biological features. This study uses publicly available FinnGen Release 9 genome-wide association study (GWAS) summary statistics and public reference resources to characterise their shared and disorder-specific common-variant genetic architecture.

The study is based entirely on summary-level and publicly available data. It does not use UK Biobank data, individual-level participant data or controlled-access biobank data.

The analytical framework includes:

1. Genome-wide genetic correlation using LD Score Regression (LDSC)
2. Local genetic correlation using LAVA
3. MAGMA gene-level association analysis
4. Exploratory pathway and gene-set annotation
5. SuSiE fine-mapping of selected shared regions
6. Single-cell RNA-seq contextualisation using published datasets
7. Cell-type heritability enrichment using stratified LDSC (S-LDSC)

## Key findings

- **Genome-wide genetic sharing:** LDSC estimated a moderate positive genetic correlation between uterine fibroids and endometriosis (rg = 0.511, SE = 0.069, P = 1.37 × 10^-13), supporting shared but incomplete common-variant architecture.
- **Regional convergence:** The primary overlap-aware LAVA analysis produced 92 valid bivariate tests and identified 16 local-correlation regions meeting the Benjamini–Hochberg false-discovery-rate threshold.
- **Gene-level architecture:** MAGMA identified 147 Bonferroni-significant genes for uterine fibroids and 38 for endometriosis. Ten genes were significant in both trait-specific analyses, while 137 were classified as fibroid-specific and 28 as endometriosis-specific under the study definition.
- **Cross-layer agreement:** Seven of the ten shared MAGMA-significant genes mapped within regions showing significant local genetic correlation.
- **Fine-mapping:** SuSiE analyses of the ESR1/SYNE1, WNT4, GREB1, WT1 and DNM3 regions supported allelic complexity at several loci, particularly ESR1/SYNE1. Fine-mapping results are interpreted as variant prioritisation rather than proof of causality.
- **Cellular context:** Single-cell contextualisation and S-LDSC provided convergent evidence for stromal involvement in both disorders, alongside disease-weighted cellular heterogeneity. These analyses do not establish a causal cell type.
- **Pathway context:** Exploratory annotations highlighted hormonal/reproductive, WNT/developmental, extracellular-matrix/fibrosis, adhesion/migration and related biological categories. These findings are hypothesis-generating.

## Data sources

### FinnGen GWAS summary statistics

- Uterine fibroids: `CD2_BENIGN_LEIOMYOMA_UTERI`
- Endometriosis: `N14_ENDOMETRIOSIS`
- Release: FinnGen R9
- Access: <https://r9.finngen.fi>

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
- GSE203191: endometriosis-related tissue

The published datasets were not re-clustered or reannotated. Original cell labels, processed expression information and ranked-marker evidence were used as described in the manuscript.

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
scripts/
  01_genetic_correlation/      LDSC and LAVA analyses
  02_gene_level/               MAGMA gene-level analysis
  03_pathway_annotation/       Exploratory pathway and gene-set annotation
  04_finemapping/              SuSiE fine-mapping
  05_single_cell_context/      Single-cell contextualisation and S-LDSC
  06_figure_generation/        Main and supplementary figure scripts

results/
  aggregate_results/           Permitted aggregate analysis outputs
  supplementary_tables/       Supplementary Tables S1-S5

figures/
  main/                        Main manuscript figures
  supplementary/              Supplementary Figures S1-S2

environment/
  R_packages.txt               R package and version information
  requirements.txt             Python package information, if applicable
```


## Suggested analysis workflow

Exact script names should match the final publication archive. A typical order is:

```bash
# 1. Genome-wide genetic correlation
Rscript scripts/01_genetic_correlation/run_ldsc_models.R

# 2. Local genetic correlation
Rscript scripts/01_genetic_correlation/prepare_lava_inputs.R
Rscript scripts/01_genetic_correlation/run_lava_blocks.R
Rscript scripts/01_genetic_correlation/process_lava_results.R

# 3. MAGMA gene-level analysis
bash scripts/02_gene_level/run_magma_gene_level.sh
# (Note: Scripts for downstream MAGMA plotting and tissue enrichment are also provided in this directory)

# 4. Exploratory pathway annotation
# (Note: Pathway annotation scripts are pending deposit in scripts/03_pathway_annotation/)

# 5. SuSiE fine-mapping
# (Note: Fine-mapping scripts are pending deposit in scripts/04_finemapping/)

# 6. Single-cell contextualisation and S-LDSC
Rscript scripts/05_single_cell_context/gse162122_fibroid_scrna_shared_gene_context.R
python scripts/05_single_cell_context/plot_scrna_shared_gene_context.py
python scripts/05_single_cell_context/combine_scrna_shared_gene_context.py

# 7. Figure generation
python scripts/06_figure_generation/generate_publication_plots.py
```


## Main manuscript outputs

### Main figures

1. Figure 1: Integrated analytical framework
2. Figure 2: Formal LAVA local genetic-correlation analysis
3. Figure 3: Gene-level architecture of uterine fibroids and endometriosis
4. Figure 4: Cell-type heritability enrichment and exploratory pathway annotation

### Main table

1. Table 1: Public data sources used in the study

### Supplementary figures

1. Supplementary Figure S1: Single-cell contextualisation of genes significant in both trait-specific MAGMA analyses
2. Supplementary Figure S2: SuSiE fine-mapping of shared genomic loci

### Supplementary tables

1. Supplementary Table S1: LDSC and LAVA analyses
2. Supplementary Table S2: Complete MAGMA gene-level results
3. Supplementary Table S3: Exploratory pathway and gene-set annotations
4. Supplementary Table S4: SuSiE fine-mapping results
5. Supplementary Table S5: Single-cell contextualisation and S-LDSC enrichment results

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
- Harmonisation, quality-control thresholds, sample-size parameters and sampling-correlation specifications should match the final manuscript and Supplementary Table S1.
- LAVA results should be reproduced using the primary LDSC-standardised overlap-aware specification. Alternative count-derived and no-overlap implementations are sensitivity analyses.
- Fine-mapping was performed separately for each disorder using ±500 kb windows, 1000 Genomes European linkage disequilibrium, a maximum of ten effects per locus and 95% credible sets.
- Pathway, fine-mapping and cellular analyses are prioritisation or contextualisation layers and should not be presented as definitive causal evidence.

## Citation

Until publication:

> Modhukur V, Lingasamy P, Patel N, Salumets A. *Integrative genomic analyses reveal regional and cellular convergence between uterine fibroids and endometriosis.* Manuscript in preparation (2026).


## Correspondence

**Vijayachitra Modhukur**  
Department of Obstetrics and Gynaecology  
Institute of Clinical Medicine, University of Tartu, Estonia  
Email: vijayachitra.modhukur@ut.ee
