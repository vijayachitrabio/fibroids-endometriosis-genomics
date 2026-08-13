# Shared genetic architecture and divergent biology of uterine fibroids and endometriosis

This repository contains analysis code, phenotype definitions, code lists, endpoint identifiers, software information and permitted aggregate/derived outputs for the manuscript:

> Modhukur V, Lingasamy P, Patel N, Salumets A. Shared genetic architecture and divergent biology of uterine fibroids and endometriosis. Manuscript in preparation (2026).

No UK Biobank individual-level data, participant identifiers, raw hospital episode records, raw diagnosis dates, raw genotype data, raw Olink matrices or other restricted participant-level files are included in this repository.

## Overview

This study integrates UK Biobank observational phenotyping with FinnGen R9 GWAS summary statistics and external molecular resources to distinguish shared susceptibility from disease-specific biology in uterine fibroids and endometriosis.
<img width="947" height="480" alt="image" src="https://github.com/user-attachments/assets/08450c43-119b-4e3e-8569-644d6b569acd" />



The analytical framework includes:

1. UK Biobank case-control and comorbidity analysis
2. Temporal ordering of recorded comorbidities relative to the index gynaecological diagnosis
3. Genome-wide genetic correlation using LDSC
4. Local genetic correlation using LAVA
5. MAGMA gene-level and exploratory pathway analysis
6. Summary-level Mendelian randomisation and exploratory PheWAS-MR
7. Bayesian colocalisation, SuSiE fine-mapping, Olink proteomics and targeted cis-pQTL protein Mendelian randomisation
8. Single-cell RNA-seq contextualisation and cell-type heritability enrichment using S-LDSC

## Key findings

- **Comorbidity:** Among 38 prespecified primary comorbidity traits, 34 were FDR-significant for uterine fibroids, 35 for endometriosis and 32 for both disorders.
- **Temporal ordering:** Among temporally evaluable traits, 9 of 16 fibroid-associated traits and 9 of 15 endometriosis-associated traits were predominantly recorded after the index diagnosis. These analyses reflect healthcare-record timing, not biological onset or causality.
- **Genome-wide genetic correlation:** LDSC estimated moderate positive genome-wide genetic correlation between uterine fibroids and endometriosis (rg = 0.511, SE = 0.069, P = 1.37 x 10^-13).
- **Local genetic correlation:** The primary LDSC-standardised overlap-aware LAVA analysis yielded 92 valid bivariate tests and 16 FDR-significant local-correlation blocks.
- **Gene-level analysis:** MAGMA identified 147 Bonferroni-significant genes for fibroids, 38 for endometriosis and 10 shared genes. Seven of the ten shared MAGMA genes mapped within significant primary LAVA regions.
- **Cellular context:** Single-cell contextualisation and S-LDSC converged on stromal biology in both disorders, while supporting disease-weighted cellular heterogeneity.
- **Proteomics:** Nineteen Olink proteins were FDR-significant for fibroids and two for endometriosis, with no shared FDR-significant protein.
- **TFPI prioritisation:** TFPI is retained only as an exploratory fibroid-weighted vascular protein-prioritisation signal. Genetically higher TFPI was associated with lower hypertension risk, but pQTL-GWAS colocalisation was non-confirmatory and mediation was not established.

## Data availability

This study used the following resources:

- **UK Biobank** — analyses conducted under UK Biobank Application 1224312. Individual-level UK Biobank data are not included in this repository and remain available only to approved researchers through UK Biobank access procedures.
- **FinnGen R9 GWAS summary statistics** — available from https://r9.finngen.fi.
- **1000 Genomes Phase 3 European reference panels** — used for linkage-disequilibrium-based analyses where appropriate.
- **Single-cell RNA-seq datasets** — fibroid/myometrium and endometriosis resources described in the manuscript and Supplementary Methods.
- **Published plasma pQTL resources** — used for targeted cis-pQTL protein Mendelian randomisation.

Only permitted aggregate, non-disclosive, publication-facing outputs are included. No participant-level UK Biobank records, participant identifiers, individual diagnosis dates, raw hospital episode data, raw genotype data or raw Olink participant-level data are redistributed.

## Repository structure

```text
scripts/
  01_comorbidity/              UK Biobank phenotype definitions, logistic regression and temporal analysis scripts
  02_genetic_correlation/      LDSC and LAVA scripts
  03_mendelian_randomisation/  Exploratory MR and PheWAS-MR scripts
  04_gene_level/               MAGMA gene-level and exploratory pathway scripts
  05_proteomics/               Olink association, cis-pQTL protein MR and pQTL-GWAS colocalisation scripts
  06_figure_generation/        Publication figure scripts
  07_single_cell_context/      Single-cell contextualisation and S-LDSC annotation scripts
  08_finemapping_coloc/        Colocalisation and SuSiE fine-mapping scripts

figures/
  Publication and supplementary figures where permitted

results/
  Permitted aggregate summary outputs only

environment/
  Python and R package information
```

## Reproducibility

The scripts are provided for transparency and reproducibility. Full reproduction of UK Biobank analyses requires approved UK Biobank access under the relevant application and access to the same controlled or restricted resources. Restricted UK Biobank participant-level data are not included and cannot be redistributed.

Before running the workflow, users must obtain the required external resources through their original access routes and comply with all relevant data-use terms.

## Suggested workflow

Run order may vary depending on local data access and file paths. A typical analysis order is:

```bash
# 1. Comorbidity and temporal analyses
Rscript scripts/01_comorbidity/precompute_results.R
Rscript scripts/01_comorbidity/run_temporal_analysis.R

# 2. Genome-wide and local genetic correlation
Rscript scripts/02_genetic_correlation/Run_LDSC_Models.R
Rscript scripts/02_genetic_correlation/run_lava_primary_overlap_aware.R

# 3. Gene-level and pathway analyses
bash scripts/04_gene_level/run_magma_gene_level.sh
Rscript scripts/04_gene_level/annotate_shared_and_specific_genes.R

# 4. Exploratory Mendelian randomisation
Rscript scripts/03_mendelian_randomisation/run_exploratory_phewas_mr.R

# 5. Colocalisation and fine-mapping
Rscript scripts/08_finemapping_coloc/run_coloc_abf.R
Rscript scripts/08_finemapping_coloc/run_susie_finemapping.R

# 6. Proteomics and protein MR
Rscript scripts/05_proteomics/run_olink_association.R
Rscript scripts/05_proteomics/run_targeted_cis_pqtl_protein_mr.R

# 7. Single-cell context and S-LDSC
Rscript scripts/07_single_cell_context/run_scrna_contextualisation.R
Rscript scripts/07_single_cell_context/run_sldsc_celltype_enrichment.R

# 8. Figure generation
Rscript scripts/06_figure_generation/generate_publication_figures.R
```

Script names may differ slightly from local development branches; the publication archive should preserve the final manuscript-specific scripts and parameters.

## Software

| Tool | Version / implementation | Reference |
|---|---|---|
| MAGMA | v1.10 | de Leeuw et al., PLoS Computational Biology 2015 |
| LDSC | LDSC / GenomicSEM workflow | Bulik-Sullivan et al., Nature Genetics 2015 |
| LAVA | v0.1.5 | Werme et al., Nature Genetics 2022 |
| coloc | coloc.abf | Giambartolomei et al., PLoS Genetics 2014 |
| susieR | SuSiE fine-mapping | Wang et al., JRSS-B 2020 |
| TwoSampleMR | MR workflow | Hemani et al., eLife 2018 |
| R | see environment/R_packages.txt |  |
| Python | see environment/requirements.txt |  |

## UK Biobank data-use note

This repository does not redistribute UK Biobank individual-level data. Any researcher wishing to reproduce UK Biobank-based analyses must apply directly to UK Biobank and comply with UK Biobank access, security, reporting and publication requirements.

Before public release, repository outputs should be checked to ensure that no restricted participant-level data, raw phenotype extracts, raw diagnosis dates, individual-level Olink values or small-cell count outputs are included.

## Citation

> Modhukur V et al. Shared genetic architecture and divergent biology of uterine fibroids and endometriosis. Manuscript in preparation (2026).

## Correspondence

Vijayachitra Modhukur — vijayachitra.modhukur@ut.ee  
Department of Obstetrics and Gynaecology, Institute of Clinical Medicine, University of Tartu, Estonia
  07_single_cell_context/ scRNA-seq cell-type validation mapping
  08_eqtl_coloc/          Bayesian colocalisation of GWAS and GTEx eQTL signals

figures/
  (Generated figures will output here)

results/
  (Result CSVs will be generated here upon running the scripts)

environment/
  requirements.txt        Python dependencies
  R_packages.txt          R package versions
```

## How to reproduce

**Prerequisites:** R ≥ 4.2, Python ≥ 3.10, MAGMA v1.10 binary (see `scripts/04_gene_level/README_MAGMA_SETUP.txt`)

Run in order:
```bash
# 1. Comorbidity analysis (R & Python)
Rscript scripts/01_comorbidity/precompute_results.R
Rscript scripts/01_comorbidity/run_temporal_analysis.R
python scripts/01_comorbidity/run_temporal_statistical_tests.py

# 2. Genetic correlation (R & Python)
Rscript scripts/02_genetic_correlation/Run_LDSC_Models.R
python scripts/02_genetic_correlation/lava_coloc_robust.py

# 3. Mendelian randomisation (R & Python)
Rscript scripts/03_mendelian_randomisation/PheWAS_MR_46traits.R
python scripts/03_mendelian_randomisation/phewas_mr_robust.py

# 4. Gene-level & Tissue Enrichment (MAGMA)
bash scripts/04_gene_level/run_tissue_enrichment.sh
python scripts/04_gene_level/plot_tissue_enrichment.py

# 5. Proteomics (R)
Rscript scripts/05_proteomics/01b_pQTL_instruments_local.R
Rscript scripts/05_proteomics/02b_coloc_local.R
Rscript scripts/05_proteomics/03b_proteome_MR_local.R
Rscript scripts/05_proteomics/04_pMR_coloc_visualization.R

# 6. eQTL Colocalization (Python)
python scripts/08_eqtl_coloc/colocalize_gwas_eqtl.py

# 7. Figures (Python & R)
python scripts/06_figure_generation/generate_publication_plots.py
python scripts/06_figure_generation/generate_global_hub_map.py
python scripts/06_figure_generation/generate_temporal_figure.py
python scripts/06_figure_generation/generate_headtohead_forest.py
Rscript scripts/03_mendelian_randomisation/generate_phewas_summary_forest.R
python scripts/06_figure_generation/generate_workflow_diagram.py

# 8. Single-cell context (R & Python)
# IMPORTANT: Before running, manually download the raw single-cell datasets 
# for Fibroids (GEO: GSE162122) and Endometriosis (GEO: GSE203191) from NCBI GEO.
Rscript scripts/07_single_cell_context/gse162122_fibroid_scrna_shared_gene_context.R
python scripts/07_single_cell_context/generate_scrna_celltype_validation_tables.py
python scripts/07_single_cell_context/plot_scrna_shared_gene_context.py

```

## Software

| Tool | Version | Reference |
|---|---|---|
| MAGMA | v1.10 | de Leeuw et al. *PLoS Comput Biol* 2015 |
| TwoSampleMR | ≥ 0.5.7 | Hemani et al. *eLife* 2018 |
| coloc | ≥ 5.1 | Giambartolomei et al. *PLoS Genet* 2014 |
| LAVA | ≥ 1.0 | Ning et al. *Nat Genet* 2021 |
| LDSC | ≥ 1.0 | Bulik-Sullivan et al. *Nat Genet* 2015 |
| Python | 3.10+ | pandas, numpy, matplotlib, scipy |
| R | 4.2+ | see environment/R_packages.txt |

## Citation

> Modhukur V et al. Shared Genetic Architecture and Temporal Comorbidity Patterns in Uterine Fibroids and Endometriosis. *Manuscript in preparation* (2026).

## Correspondence

Vijayachitra Modhukur — vijayachitra.modhukur@celvia.ee  
Institute of Clinical Medicine, University of Tartu, Estonia | Celvia CC, Tartu, Estonia
