#!/usr/bin/env Rscript

# This script generates the standardized LDSC overlap matrix for LAVA.
# The genetic covariance intercept is 0.1917, and the LDSC intercepts are 1.0967 and 1.1274.
# The sampling correlation is calculated as: covariance_intercept / sqrt(intercept1 * intercept2)

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x
args_full <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE)[1] %||% "")
root <- if (nzchar(self)) {
  normalizePath(file.path(dirname(self), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
if (!dir.exists(file.path(root, "outputs_v2"))) root <- normalizePath(".", mustWork = FALSE)

out_dir <- file.path(root, "outputs_v2", "formal_ld_audit", "lava")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# LDSC intercepts and covariance intercept
cov_intercept <- 0.1917
int1 <- 1.0967
int2 <- 1.1274

# Calculate sampling correlation
sampling_cor <- cov_intercept / sqrt(int1 * int2)

# Create matrix (LAVA expects it in this format)
overlap_matrix <- matrix(
  c(1, sampling_cor, sampling_cor, 1),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Endometriosis", "Uterine_Fibroids"), c("Endometriosis", "Uterine_Fibroids"))
)

overlap_file <- file.path(out_dir, "lava_sample_overlap_ldsc_standardized.tsv")
write.table(overlap_matrix, overlap_file, sep = "\t", quote = FALSE, col.names = NA)
message("Saved standardized LDSC overlap matrix to: ", overlap_file)
