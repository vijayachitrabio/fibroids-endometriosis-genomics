#!/usr/bin/env Rscript
# ============================================================
# Process formal LAVA 2,495-block overlap-specified results → replacement S2 table
# Run after run_lava_full_2495.sh completes:
#   Rscript formal_ld_workflows/process_lava_results.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

lava_dir <- "outputs_v2/formal_ld_audit/lava"

# Use the overlap-specified formal LAVA run as primary. The no-overlap run is
# retained only as an archived sensitivity analysis.
results_file <- file.path(lava_dir, "LAVA_Formal_2495_Blocks_Results_overlap_proxy_lava_loci_2495_blocks.csv")
stopifnot(file.exists(results_file))

cat("Reading results:", results_file, "\n")
res <- fread(results_file)
cat("Total blocks with bivariate result:", nrow(res), "\n")

# Remove rows with missing p-value
res <- res[!is.na(p)]
cat("Blocks with valid p-value:", nrow(res), "\n")

# FDR correction
res[, p_fdr := p.adjust(p, method = "fdr")]

# Significant at FDR < 0.05
sig <- res[p_fdr < 0.05]
cat("FDR < 0.05 blocks:", nrow(sig), "\n")

# Sort by p-value
sig <- sig[order(p)]

# Round for readability
round_cols <- c("rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper")
for (col in intersect(round_cols, names(sig))) {
  sig[[col]] <- round(sig[[col]], 4)
}
sig[, p := signif(p, 3)]
sig[, p_fdr := signif(p_fdr, 3)]

# Save S2 replacement table
out_s2 <- file.path(lava_dir, "Table_S2_LAVA_Local_rg_Formal_2495Blocks_overlap_proxy_FDR05.csv")
fwrite(sig, out_s2)
cat("Saved S2 replacement table:", out_s2, "\n")
cat("  Rows:", nrow(sig), "\n")
cat("\nTop 10 significant blocks:\n")
print(sig[1:min(10, nrow(sig)), .(LOC, CHR, START, STOP, rho, p, p_fdr)])

# Also save full results with FDR column
out_full <- file.path(lava_dir, "Table_LAVA_Full_2495Blocks_AllResults_overlap_proxy_withFDR.csv")
fwrite(res[order(p)], out_full)
cat("\nSaved full results with FDR:", out_full, "\n")

# Copy to manuscript assets
manu_dir <- "manuscript_assets_latest_2026_05_20/supplementary_tables"
if (dir.exists(manu_dir)) {
  file.copy(out_s2, file.path(manu_dir, "Table_S2_LAVA_Local_rg_Formal_2495Blocks_overlap_proxy_FDR05.csv"), overwrite = TRUE)
  cat("Copied to manuscript assets:", manu_dir, "\n")
}

cat("\n=== Summary ===\n")
cat("Blocks with valid p-values:", nrow(res), "\n")
cat("Blocks with FDR < 0.05:", nrow(sig), "\n")
if (nrow(sig) > 0) {
  cat("rho range in significant blocks:", round(min(sig$rho, na.rm=TRUE), 3),
      "to", round(max(sig$rho, na.rm=TRUE), 3), "\n")
}
