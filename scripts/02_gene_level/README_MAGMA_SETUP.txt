============================================================
  MAGMA Gene-Level Analysis — Setup & Run Instructions
  Uterine Fibroids and Endometriosis GWAS
  Prepared: April 2026
============================================================

ALL INPUTS ARE READY. You only need to:
  (1) Download the MAGMA binary
  (2) Copy two files to this folder
  (3) Run one command

------------------------------------------------------------
WHAT IS ALREADY PREPARED IN THIS FOLDER
------------------------------------------------------------

fibroids.snploc       12,018,721 SNPs with CHR/BP coordinates
                       (merged from Cleaned_Uterine_Fibroids.txt
                        × g1000_eur.bim; 63.9% match rate)

fibroids.pval         Same SNPs with P-values and N=218,728

endometriosis.snploc  12,015,108 SNPs with CHR/BP coordinates
endometriosis.pval    Same SNPs with P-values and N=218,728

run_magma.sh          Complete MAGMA pipeline (annotate + gene test)
plot_magma_results.py Post-processing: figures, tables, manuscript tables

------------------------------------------------------------
STEP 1 — DOWNLOAD MAGMA BINARY (do this on your local machine)
------------------------------------------------------------

Go to:  https://ctg.cncr.nl/software/magma

Download the Linux static binary:
  magma_v1.10_static.zip  (or latest version)

Unzip to get:
  magma                 (the executable)
  NCBI37.3.gene.loc     (gene location file, GRCh37/hg19)
  aux_files/            (optional)

NOTE: If on macOS, download the macOS version instead.
      The pipeline script works for both.

------------------------------------------------------------
STEP 2 — COPY TWO FILES TO THIS FOLDER
------------------------------------------------------------

Copy these two files into:
  uterine_fibroids/magma_inputs/

  → magma               (the MAGMA executable)
  → NCBI37.3.gene.loc   (the gene location file)

After copying, this folder should contain:
  magma
  NCBI37.3.gene.loc
  fibroids.snploc
  fibroids.pval
  endometriosis.snploc
  endometriosis.pval
  run_magma.sh
  plot_magma_results.py
  README_MAGMA_SETUP.txt

------------------------------------------------------------
STEP 3 — RUN THE ANALYSIS
------------------------------------------------------------

Open a terminal, navigate to this folder, and run:

  chmod +x magma
  bash run_magma.sh

The script will:
  1. Check all required files are present
  2. Annotate SNPs to genes (window 0 kb)
  3. Run gene-level test using 1000G EUR LD reference
     (located at: ../mr_analysis_2026_04_16/reference_data/g1000_eur)
  4. Print a summary of significant genes

Runtime:  approximately 45–90 minutes per trait

------------------------------------------------------------
STEP 4 — GENERATE FIGURES AND TABLES
------------------------------------------------------------

After run_magma.sh completes:

  python3 plot_magma_results.py

This produces:
  results/fibroids_genes.genes.out
  results/endometriosis_genes.genes.out
  results/fibroids_vs_endo_comparison.tsv

  outputs_v2/Figures_Publication/Figure4_MAGMA_GeneLevel.png
  outputs_v2/Figures_Publication/FigS7_MAGMA_Gene_Level.png
  outputs_v2/Table_MAGMA_TopGenes.csv
  outputs_v2/Supplementary_Table_S14a_Fibroids_MAGMA_AllGenes.tsv
  outputs_v2/Supplementary_Table_S14b_Endometriosis_MAGMA_AllGenes.tsv
  outputs_v2/Supplementary_Table_S14c_Shared_Genes_MAGMA.tsv

------------------------------------------------------------
TECHNICAL NOTES
------------------------------------------------------------

LD Reference:  1000 Genomes Phase 3 EUR (g1000_eur)
               22,665,064 SNPs, 503 samples
               Located at: ../mr_analysis_2026_04_16/reference_data/

Genome build:  GRCh37 / hg19 (matches NCBI37.3.gene.loc)

SNP match rate explanation:
  ~12M out of 18.8M GWAS SNPs matched to the 1000G bim file.
  The ~6.8M unmatched SNPs are typically:
    - Rare variants (MAF<1% in 1000G EUR)
    - Indels not in 1000G
    - Non-rs SNPs
  This is normal and does not affect the analysis — MAGMA
  uses the matched SNPs for gene-level aggregation.

MAGMA method (mean model):
  Aggregates SNP-level chi-squared statistics within each gene
  using a permutation-based approach that correctly accounts
  for LD. This is superior to the Fisher's combination method
  previously used in the manuscript (which does NOT correct
  for LD and is anti-conservative).

Bonferroni threshold:  P < 2.5×10⁻⁶  (0.05 / 19,993 genes)

------------------------------------------------------------
QUESTIONS?
------------------------------------------------------------
If the magma binary fails to execute, check:
  file magma           → should say "ELF 64-bit" or similar
  uname -m             → should be x86_64 (for Linux binary)
  ./magma --help       → should print MAGMA help text
============================================================
