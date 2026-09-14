#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(LAVA)
})

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

args_full <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE)[1] %||% "")
root <- if (nzchar(self)) normalizePath(file.path(dirname(self), ".."), mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
if (!dir.exists(file.path(root, "gwas_scripts"))) root <- normalizePath(".", mustWork = TRUE)

ref_prefix <- file.path(root, "mr_analysis_2026_04_16", "reference_data", "g1000_eur")
in_dir <- file.path(root, "outputs_v2", "formal_ld_audit", "lava")
out_file <- file.path(in_dir, "LAVA_Formal_2495_Blocks_Results.csv")
progress_file <- file.path(in_dir, "LAVA_Formal_2495_Blocks_Progress.csv")
qc_file <- file.path(in_dir, "LAVA_Formal_2495_Blocks_QC_Log.csv")

input_info <- file.path(in_dir, "lava_input_info.tsv")
overlap_file <- file.path(in_dir, "lava_sample_overlap_ldsc_corr.tsv")
loci_file <- get_arg("loci-file", file.path(in_dir, "lava_loci_2495_blocks.tsv"))
stopifnot(file.exists(input_info), file.exists(overlap_file), file.exists(loci_file))
stopifnot(file.exists(paste0(ref_prefix, ".bed")), file.exists(paste0(ref_prefix, ".bim")), file.exists(paste0(ref_prefix, ".fam")))

max_loci <- as.integer(get_arg("max-loci", NA))
start_at <- as.integer(get_arg("start-at", 1))
use_overlap <- get_arg("overlap", "yes") != "no"
sample_overlap <- if (use_overlap) overlap_file else NULL
mode_suffix <- if (use_overlap) "overlap_ldsc_corr" else "no_overlap"
mode_suffix <- paste0(mode_suffix, "_", tools::file_path_sans_ext(basename(loci_file)))

out_file <- file.path(in_dir, paste0("LAVA_Formal_2495_Blocks_Results_", mode_suffix, ".csv"))
progress_file <- file.path(in_dir, paste0("LAVA_Formal_2495_Blocks_Progress_", mode_suffix, ".csv"))
qc_file <- file.path(in_dir, paste0("LAVA_Formal_2495_Blocks_QC_Log_", mode_suffix, ".csv"))

message("Processing LAVA input...")
input <- process.input(
  input.info.file = input_info,
  sample.overlap.file = sample_overlap,
  ref.prefix = ref_prefix,
  input.dir = in_dir
)

loci <- read.loci(loci_file)
if (!is.na(max_loci)) loci <- loci[seq_len(min(max_loci, nrow(loci))), ]
if (start_at > 1) loci <- loci[seq(from = start_at, to = nrow(loci)), ]

done <- character()
if (file.exists(progress_file)) {
  existing <- fread(progress_file)
  if ("LOC" %in% names(existing)) done <- unique(existing$LOC)
}

append_row <- function(path, row) {
  fwrite(row, path, append = file.exists(path))
}

message("Running formal LAVA over ", nrow(loci), " blocks...")
for (i in seq_len(nrow(loci))) {
  loc <- loci[i, ]
  loc_id <- as.character(loc$LOC)
  if (loc_id %in% done) next
  if (i %% 25 == 0) message("  processed attempt ", i, " / ", nrow(loci))

  qc_base <- data.table(
    LOC = loc_id,
    CHR = loc$CHR,
    START = loc$START,
    STOP = loc$STOP
  )

  locus <- tryCatch(
    process.locus(loc, input, drop.failed = FALSE, max.block.size = 3000),
    error = function(e) {
      append_row(qc_file, cbind(qc_base, status = "process_locus_error", message = conditionMessage(e)))
      NULL
    }
  )
  if (is.null(locus)) next

  res <- tryCatch(
    run.univ.bivar(locus, univ.thresh = 1, p.values = TRUE, CIs = TRUE),
    error = function(e) {
      append_row(qc_file, cbind(qc_base, status = "run_univ_bivar_error", message = conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res) || is.null(res$biv) || nrow(res$biv) == 0) {
    append_row(qc_file, cbind(qc_base, status = "no_bivariate_result", message = "No bivariate result returned"))
    next
  }

  biv <- as.data.table(res$biv)
  univ <- if (!is.null(res$univ)) as.data.table(res$univ) else data.table()
  out <- cbind(
    qc_base,
    n_snps = locus$n.snps,
    K = locus$K,
    biv
  )
  if (nrow(univ) > 0) {
    # Keep the raw univariate table in QC rather than forcing fragile wide names.
    append_row(qc_file, cbind(qc_base, status = "univ", message = paste(capture.output(print(univ)), collapse = " | ")))
  }
  append_row(progress_file, out)
}

if (file.exists(progress_file)) {
  final <- fread(progress_file)
  fwrite(final, out_file)
  message("Saved formal LAVA results: ", out_file)
} else {
  message("No formal LAVA results were produced. Check QC log: ", qc_file)
}
