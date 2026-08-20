#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

args_full <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE)[1] %||% "")
root <- if (nzchar(self)) {
  normalizePath(file.path(dirname(self), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
if (!dir.exists(file.path(root, "gwas_scripts"))) root <- normalizePath(".", mustWork = TRUE)

gwas_dir <- file.path(root, "gwas_scripts")
ref_dir <- file.path(root, "mr_analysis_2026_04_16", "reference_data")
out_dir <- file.path(root, "outputs_v2", "formal_ld_audit", "lava")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- data.table(
  phenotype = c("Endometriosis", "Uterine_Fibroids"),
  file = file.path(
    gwas_dir,
    c("finngen_R9_N14_ENDOMETRIOSIS.gz", "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz")
  ),
  cases = c(17045L, 35474L),
  controls = c(239335L - 17045L, 218728L - 35474L),
  n_total = c(239335L, 218728L)
)

ref_bim <- file.path(ref_dir, "g1000_eur.bim")
blocks_in <- file.path(ref_dir, "LAVA_s2500_m25_f1_w200.blocks")
stopifnot(file.exists(ref_bim), file.exists(blocks_in))
stopifnot(all(file.exists(traits$file)))

message("Reading 1000G EUR reference SNP IDs...")
ref <- fread(ref_bim, select = c(2, 5, 6), col.names = c("SNP", "REF_A1", "REF_A2"))
ref[, SNP := tolower(SNP)]
ref_ids <- ref$SNP

prepare_trait <- function(row) {
  pheno <- row$phenotype
  out <- file.path(out_dir, paste0(pheno, "_lava_sumstats.tsv.gz"))
  if (file.exists(out)) {
    message("Already prepared: ", out)
    return(out)
  }

  message("Preparing LAVA summary statistics for ", pheno, "...")
  dt <- fread(
    row$file,
    select = c("rsids", "alt", "ref", "beta", "sebeta", "pval"),
    showProgress = TRUE
  )
  setnames(dt, c("rsids", "alt", "ref", "beta", "sebeta", "pval"),
           c("SNP_RAW", "A1", "A2", "BETA", "SE", "P"))
  dt[, SNP := tolower(tstrsplit(SNP_RAW, ",", fixed = TRUE, keep = 1L)[[1]])]
  dt <- dt[SNP != "" & !is.na(SNP)]
  dt <- dt[SNP %chin% ref_ids]
  dt <- dt[!duplicated(SNP)]
  dt <- dt[
    !is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(P) & P > 0 & P <= 1 &
      nchar(A1) == 1 & nchar(A2) == 1,
    .(SNP, A1, A2, BETA, P, N = row$n_total)
  ]
  fwrite(dt, out, sep = "\t")
  message("  wrote ", nrow(dt), " SNPs: ", out)
  out
}

sumstat_files <- vapply(seq_len(nrow(traits)), function(i) prepare_trait(traits[i]), character(1))

input_info <- traits[, .(
  phenotype,
  cases,
  controls,
  filename = basename(sumstat_files)
)]
input_info_file <- file.path(out_dir, "lava_input_info.tsv")
fwrite(input_info, input_info_file, sep = "\t")

# FinnGen phenotypes share participants/controls. Use the LDSC-derived overlap
# proxy from the audit table as the primary formal setting, and run a NULL-overlap
# sensitivity in the LAVA script if desired.
overlap_n <- 59264
overlap <- matrix(
  c(traits$n_total[1], overlap_n, overlap_n, traits$n_total[2]),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(traits$phenotype, traits$phenotype)
)
overlap_file <- file.path(out_dir, "lava_sample_overlap_counts.tsv")
write.table(overlap, overlap_file, sep = "\t", quote = FALSE, col.names = NA)

blocks <- fread(blocks_in)
setnames(blocks, tolower(names(blocks)))
blocks[, LOC := sprintf("block_%04d", seq_len(.N))]
loci <- blocks[, .(LOC, CHR = chr, START = start, STOP = stop)]
loci_file <- file.path(out_dir, "lava_loci_2495_blocks.tsv")
fwrite(loci, loci_file, sep = "\t")

message("Prepared formal LAVA inputs:")
message("  ", input_info_file)
message("  ", overlap_file)
message("  ", loci_file)
