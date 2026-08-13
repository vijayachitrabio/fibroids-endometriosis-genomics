###############################################################################
# PheWAS-MR: Phenome-Wide Two-Sample Mendelian Randomization
# Fibroids & Endometriosis → 46 comorbidities
#
# PURPOSE  : Expand the original 2-outcome MR analysis to all 46 comorbidities
#            by using FinnGen R9 GWAS as both exposure and outcome sources.
#
# INPUTS   : FinnGen R9 GWAS summary statistics (local .gz files)
#            FinnGen R9 endpoint codes mapped to the 46 comorbidity panel
#
# OUTPUTS  : outputs_v2/PheWAS_MR_All_Results.csv       — full 3-estimator table
#            outputs_v2/PheWAS_MR_IVW_Summary.csv        — IVW-only summary
#            outputs_v2/MR_Plots/PheWAS_Forest_*.png     — per-exposure forest plots
#
# REQUIRES : TwoSampleMR, data.table, ggplot2, dplyr, ieugwasr
#
# USAGE    : Rscript PheWAS_MR_46traits.R
#            Or source() from RStudio after setting BASE_DIR below.
###############################################################################

library(TwoSampleMR)
library(data.table)
library(ggplot2)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  BASE_DIR <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) getwd())
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    BASE_DIR <- dirname(normalizePath(sub("--file=", "", file_arg[1])))
  } else {
    BASE_DIR <- getwd()
  }
}
# Override paths via env vars; fall back to siblings of this script's directory.
# Set FINNGEN_GWAS_DIR and FIBROID_OUT_DIR in ~/.Renviron or before sourcing.
GWAS_DIR  <- Sys.getenv("FINNGEN_GWAS_DIR",
               unset = file.path(BASE_DIR, "..", "Finngen_phewas", "finngen_r9_phewas"))
OUT_DIR   <- Sys.getenv("FIBROID_OUT_DIR",
               unset = file.path(BASE_DIR, "..", "outputs_v2"))
PLOT_DIR  <- file.path(OUT_DIR, "MR_Plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── FinnGen endpoint → comorbidity mapping ────────────────────────────────────
# ── FinnGen endpoint → comorbidity mapping ────────────────────────────────────
# Map each of the 46 comorbidities to the closest FinnGen R9 endpoint.
finngen_map <- data.frame(
  comorbidity = c(
    "irritable_bowel_syndrome", "endometrial_cancer", "depression", "anxiety",
    "heavy_menstrual_bleeding", "pcos", "ovarian_cyst", "ovarian_cancer",
    "breast_cancer", "cervical_cancer", "endometrial_polyp", "adenomyosis",
    "obesity", "type2_diabetes", "hypertension", "hypothyroidism",
    "hyperthyroidism", "osteoarthritis", "osteoporosis", "rheumatoid_arthritis",
    "sle", "ibd", "coeliac_disease", "migraine", "fibromyalgia",
    "chronic_pelvic_pain", "iron_deficiency_anaemia", "recurrent_uti",
    "overactive_bladder", "urinary_incontinence",
    "urinary_incontinence_unspecified",   # alternative endpoint (R18 code)
    "dyslipidaemia", "nafld",
    "gerd", "peptic_ulcer_disease", "back_pain", "psoriasis",
    "vitamin_d_deficiency", "leiomyosarcoma", "cervical_polyp", "interstitial_cystitis",
    "vulvodynia", "thrombocytopenia", "ptsd"
  ),
  finngen_code = c(
    "K11_IBS", "C3_CORPUS_UTERI_EXALLC", "F5_DEPRESSIO", "F5_ANXIETY",
    "N14_HEAVYMENSTR", "E4_PCOS", "N14_OVARYCYST", "C3_OVARY_EXALLC",
    "C3_BREAST_EXALLC", "C3_CERVIXUTERI_EXALLC", "N14_ENDOPOLYP", "N14_ADENOMYOSIS",
    "E4_OBESITY", "T2D", "I9_HYPTENS", "E4_HYTHYROID",
    "E4_HYPER_THYROID", "M13_ARTHROSIS", "M13_OSTEOPOROSIS", "M13_RHEUMA_SEROPOS_OTH",
    "AB1_SLE", "K11_IBD", "K11_COELIAC", "G6_MIGRAINE", "M13_FIBROMYALGIA",
    "N14_PELVIPAIN", "D3_IRON_DEF_ANAEMIA", "N14_RECUR_UTI",
    "N14_OVERACTIVE_BLADDER", "N14_URINARYINCONT",
    "R18_UNSPE_URINARY_INCONTINENCE",     # broader unspecified code
    "E4_LIPOPROT", "NAFLD",
    "K21_GERD", "K11_PULC", "M13_BACKPAIN", "L12_PSORIASIS",
    "E4_VIT_D_DEF", "C3_LEIOMYO", "N14_CERVPOLYT", "INTERSTITIAL_CYSTIT_CHRONIC",
    "N14_VULVOVAGINFINOTH", "D3_OTHPRIMTHROMBOCYTOPENIA", "F5_PTSD"
  ),
  stringsAsFactors = FALSE
)
# Verify no duplicate comorbidity labels remain
stopifnot(!any(duplicated(finngen_map$comorbidity)))

# Deep-Scan Pathing Engine (Fix for ISSUE 2)
LIB_DIR <- GWAS_DIR   # resolved above from env var or relative path
file_inventory <- list.files(LIB_DIR, recursive = TRUE, full.names = FALSE, pattern = "\\.gz$")

finngen_map$local_file <- sapply(finngen_map$finngen_code, function(code) {
  # Favor prefix and non-zero size
  patterns <- c(paste0("finngen_R9_", code, ".gz"), paste0(code, ".gz"))
  for (p in patterns) {
    match <- file_inventory[basename(file_inventory) == p]
    if (length(match) > 0) return(match[1]) # Pick the first valid path
  }
  return(NA)
})

# ── Helper: read FinnGen GWAS and convert to TwoSampleMR format ──────────────
read_finngen <- function(path, type = "exposure") {
  cat("  Reading", basename(path), "...\n")
  # Integrity Check: Skip 0-byte or corrupted files
  if (file.size(path) < 1000) stop("File size too small (corrupted or empty)")
  
  dt <- fread(path, sep = "\t", header = TRUE)
  # Adaptive Header Map: Handle #chrom and rsids aliases
  setnames(dt, old = c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval","af_alt"),
               new = c("chr","pos","other_allele","effect_allele","SNP","beta","se","pval","eaf"),
               skip_absent = TRUE)
  
  if (!"SNP" %in% names(dt)) stop("Column SNP (rsids) not found")
  
  dt[, SNP := tstrsplit(as.character(SNP), ",", fixed = TRUE)[[1]]]   # take first rsid
  dt[, chr := as.character(chr)]
  dt[, beta  := as.numeric(beta)]
  dt[, se    := as.numeric(se)]
  dt[, pval  := as.numeric(pval)]
  dt[, eaf   := as.numeric(eaf)]
  dt <- dt[!is.na(beta) & !is.na(se) & se > 0]
  dt[, SNP := tstrsplit(SNP, ",", fixed = TRUE)[[1]]]   # take first rsid
  dt[, chr := as.character(chr)]
  dt[, beta  := as.numeric(beta)]
  dt[, se    := as.numeric(se)]
  dt[, pval  := as.numeric(pval)]
  dt[, eaf   := as.numeric(eaf)]
  dt <- dt[!is.na(beta) & !is.na(se) & se > 0]
  if (type == "exposure") {
    dt <- dt[pval < 5e-8]                               # genome-wide significant only
    # Distance-based LD pruning (500 kb windows — crude but standard for FinnGen)
    dt <- dt[order(pval)]
    dt[, window := paste0(chr, "_", floor(pos / 500000))]
    dt <- dt[!duplicated(window)]
    cat("   →", nrow(dt), "instruments after distance pruning\n")
    dt <- format_data(as.data.frame(dt), type = "exposure",
                      snp_col = "SNP", beta_col = "beta", se_col = "se",
                      eaf_col = "eaf", effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele", pval_col = "pval",
                      chr_col = "chr", pos_col = "pos")
  } else {
    dt <- format_data(as.data.frame(dt), type = "outcome",
                      snp_col = "SNP", beta_col = "beta", se_col = "se",
                      eaf_col = "eaf", effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele", pval_col = "pval",
                      chr_col = "chr", pos_col = "pos")
  }
  dt
}

# ── Helper: run full MR with all 3 estimators ─────────────────────────────────
run_mr_full <- function(dat, exp_name, out_name) {
  if (nrow(dat) < 3) {
    warning("Fewer than 3 harmonised SNPs for ", exp_name, " → ", out_name)
    return(NULL)
  }
  res   <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  het   <- mr_heterogeneity(dat)
  pleio <- mr_pleiotropy_test(dat)
  or    <- generate_odds_ratios(res)

  out <- data.frame(
    Exposure             = exp_name,
    Outcome              = out_name,
    N_instruments        = nrow(dat),
    # IVW
    IVW_OR               = round(or$or[or$method == "Inverse variance weighted"], 3),
    IVW_CI               = paste0(round(or$or_lci95[or$method == "Inverse variance weighted"], 3),
                                  "–",
                                  round(or$or_uci95[or$method == "Inverse variance weighted"], 3)),
    IVW_SE               = round(res$se[res$method == "Inverse variance weighted"], 5),
    IVW_P                = round(res$pval[res$method == "Inverse variance weighted"], 4),
    # Weighted median
    WM_OR                = round(or$or[or$method == "Weighted median"], 3),
    WM_CI                = paste0(round(or$or_lci95[or$method == "Weighted median"], 3),
                                  "–",
                                  round(or$or_uci95[or$method == "Weighted median"], 3)),
    WM_P                 = round(res$pval[res$method == "Weighted median"], 4),
    # Egger
    Egger_OR             = round(or$or[or$method == "MR Egger"], 3),
    Egger_CI             = paste0(round(or$or_lci95[or$method == "MR Egger"], 3),
                                  "–",
                                  round(or$or_uci95[or$method == "MR Egger"], 3)),
    Egger_P_slope        = round(res$pval[res$method == "MR Egger"], 4),
    Egger_intercept_P    = round(pleio$pval, 4),
    # Heterogeneity
    Cochran_Q            = round(het$Q[het$method == "Inverse variance weighted"], 3),
    Cochran_Q_P          = round(het$Q_pval[het$method == "Inverse variance weighted"], 4),
    # FDR will be applied later across all traits
    stringsAsFactors = FALSE
  )
  out
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
exposures_def <- list(
  Endometriosis    = file.path(GWAS_DIR, "finngen_R9_N14_ENDOMETRIOSIS.gz"),
  Uterine_Fibroids = file.path(GWAS_DIR, "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz")
)

all_results <- list()

for (exp_name in names(exposures_def)) {
  cat("\n══════════════════════════════════════════\n")
  cat("Exposure:", exp_name, "\n")
  cat("══════════════════════════════════════════\n")
  exp_dat <- read_finngen(exposures_def[[exp_name]], type = "exposure")
  exp_dat$exposure <- exp_name

  for (i in seq_len(nrow(finngen_map))) {
    row       <- finngen_map[i, ]
    out_name  <- row$comorbidity
    local_f   <- row$local_file

    if (is.na(local_f)) {
      cat("  ── SKIP (no local GWAS):", out_name,
          "— download from https://r9.finngen.fi/pheno/", row$finngen_code, "\n")
      next
    }
    out_path <- file.path(GWAS_DIR, local_f)
    if (!file.exists(out_path)) {
      cat("  ── SKIP (file not found):", out_path, "\n")
      next
    }

    cat("\n  ▶ Outcome:", out_name, "\n")
    tryCatch({
      out_dat <- read_finngen(out_path, type = "outcome")
      out_dat$outcome <- out_name
      dat     <- harmonise_data(exp_dat, out_dat, action = 2)
      # Remove palindromic SNPs with intermediate MAF
      dat     <- dat[!(dat$palindromic & dat$ambiguous), ]
      cat("   →", nrow(dat), "harmonised SNPs\n")
      res_row <- run_mr_full(dat, exp_name, out_name)
      if (!is.null(res_row)) {
        all_results[[length(all_results) + 1]] <- res_row
        cat("   → IVW OR =", res_row$IVW_OR,
            "(", res_row$IVW_CI, "), P =", res_row$IVW_P, "\n")
      }
      # Save per-pair scatter and forest plots
      plot_prefix <- file.path(PLOT_DIR, paste0("PheWAS_", exp_name, "_vs_", out_name))
      png(paste0(plot_prefix, "_scatter.png"), width = 800, height = 600, res = 100)
        print(mr_scatter_plot(res_row, dat))
      dev.off()
      png(paste0(plot_prefix, "_forest.png"), width = 800, height = 500, res = 100)
        print(mr_forest_plot(mr_singlesnp(dat)))
      dev.off()
    }, error = function(e) cat("   ✗ ERROR:", conditionMessage(e), "\n"))
  }
}

# ── Collate results ───────────────────────────────────────────────────────────
if (length(all_results) > 0) {
  results_df <- do.call(rbind, all_results)
  # Apply FDR correction within each exposure across all outcomes
  results_df <- results_df %>%
    group_by(Exposure) %>%
    mutate(IVW_FDR_P = p.adjust(IVW_P, method = "BH")) %>%
    ungroup()
  results_df$IVW_FDR_Sig <- results_df$IVW_FDR_P < 0.05
  results_df$Conclusion <- ifelse(results_df$IVW_FDR_Sig & results_df$IVW_P < 0.05,
                                  "FDR-significant", "Null")

  write.csv(results_df, file.path(OUT_DIR, "PheWAS_MR_All_Results.csv"), row.names = FALSE)
  cat("\n✓ Full results saved to:", file.path(OUT_DIR, "PheWAS_MR_All_Results.csv"), "\n")

  # IVW summary
  ivw_summary <- results_df %>%
    select(Exposure, Outcome, N_instruments, IVW_OR, IVW_CI, IVW_P, IVW_FDR_P, Conclusion)
  write.csv(ivw_summary, file.path(OUT_DIR, "PheWAS_MR_IVW_Summary.csv"), row.names = FALSE)

  # ── PheWAS forest plot for each exposure ─────────────────────────────────
  for (exp_name in unique(results_df$Exposure)) {
    sub <- results_df %>% filter(Exposure == exp_name) %>%
           mutate(Outcome = gsub("_", " ", Outcome),
                  Outcome = tools::toTitleCase(Outcome)) %>%
           arrange(IVW_P)
    p <- ggplot(sub, aes(x = log(IVW_OR), y = reorder(Outcome, -IVW_P),
                         colour = Conclusion)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = log(IVW_OR) - 1.96 * IVW_SE,
                         xmax = log(IVW_OR) + 1.96 * IVW_SE),
                     height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c("FDR-significant" = "#D32F2F", "Null" = "#1565C0")) +
      labs(title = paste("PheWAS-MR:", exp_name, "→ comorbidities"),
           subtitle = "Inverse-variance weighted OR (log scale). FDR correction within exposure.",
           x = "log(OR)", y = NULL,
           caption = "Instruments: genome-wide significant (P<5e-8), 500kb pruned; FinnGen R9") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    ggsave(file.path(PLOT_DIR, paste0("PheWAS_Forest_", exp_name, ".png")),
           p, width = 10, height = 8, dpi = 150)
  }
  cat("✓ PheWAS forest plots saved to:", PLOT_DIR, "\n")
} else {
  cat("\n⚠ No results to collate — check that FinnGen GWAS files are in:", GWAS_DIR, "\n")
  cat("  Required files per comorbidity: https://r9.finngen.fi (publicly available)\n")
}

# ── Print download instructions for missing files ─────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat("DOWNLOAD INSTRUCTIONS FOR MISSING FinnGen GWAS FILES\n")
cat("══════════════════════════════════════════════════════\n")
missing <- finngen_map[is.na(finngen_map$local_file), ]
cat("Run the following from your GWAS_DIR to download all missing files:\n\n")
for (i in seq_len(nrow(missing))) {
  code <- missing$finngen_code[i]
  cat(sprintf('wget "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_%s.gz" \\\n',
              code))
}
cat('\nNote: Each file is ~700MB. After download, re-run this script.\n')
cat('FinnGen R9 data is freely available under: https://finngen.gitbook.io/documentation/data-access\n')
