# MR_phenotype_divergence.R
# Two-Sample Mendelian Randomization Script for Fibroids & Endometriosis

# ==============================================================================
# OPENGWAS API TOKEN — loaded securely from environment (never hardcode in script)
#
# SETUP (one-time):
#   1. Get / renew your token at: https://api.opengwas.io/profile/
#   2. Add this line to ~/.Renviron  (RStudio: usethis::edit_r_environ())
#        OPENGWAS_JWT=eyJhbGci...your_token_here...
#   3. Restart R / RStudio
#
# SECURITY: If a token was ever committed to git, revoke it immediately at
#   https://api.opengwas.io/account before pushing any new commits.
# ==============================================================================
Sys.setenv(OPENGWAS_JWT = Sys.getenv("OPENGWAS_JWT"))
if (nchar(Sys.getenv("OPENGWAS_JWT")) == 0)
  stop("OPENGWAS_JWT is not set. Add it to ~/.Renviron and restart R.")
# ==============================================================================

library(ieugwasr)
# We wrap this in a suppress warning to deal with environments where it might load slowly
suppressWarnings(suppressPackageStartupMessages(library(TwoSampleMR)))
library(data.table)

out_folder <- "../outputs_v2"
dir.create(out_folder, showWarnings = FALSE)

cat("Starting Mendelian Randomization Pipeline...\n")

# Define Exposures (FinnGen R9 equivalents/stable IEU IDs)
exposure_ids <- list(
  Endometriosis    = "finn-b-N14_ENDOMETRIOSIS",
  Uterine_Fibroids = "finn-b-D3_CODE_FIBROMA_UTERI"
)

# Define Outcomes
outcome_ids <- list(
  IBS                    = "finn-b-K11_IBS",
  Heavy_Menstrual_Bleeding= "finn-b-N14_EXCESMENSTR",
  Endometrial_Cancer     = "finn-b-C3_ENDOMETRIUM",
  Type2_Diabetes         = "ebi-a-GCST006867"
)

all_results <- list()

for (exp_name in names(exposure_ids)) {
  exp_id <- exposure_ids[[exp_name]]
  cat("\n=======================================================\n")
  cat("Extracting Instruments for:", exp_name, "(", exp_id, ")...\n")
  
  # Step 1: Extract Instruments (SNPs strongly associated with exposure)
  exp_data <- extract_instruments(outcomes = exp_id)
  
  if (is.null(exp_data) || nrow(exp_data) == 0) {
    cat("Warning: No instruments found for", exp_name, "Skipping...\n")
    next
  }
  
  exp_data$exposure <- exp_name
  cat("Found", nrow(exp_data), "SNPs.\n")
  
  for (out_name in names(outcome_ids)) {
    out_id <- outcome_ids[[out_name]]
    cat("\n  Testing Outcome:", out_name, "(", out_id, ")...\n")
    
    # Step 2: Extract Outcome Data
    out_data <- extract_outcome_data(
      snps = exp_data$SNP,
      outcomes = out_id,
      proxies = FALSE
    )
    
    if (is.null(out_data) || nrow(out_data) == 0) {
      cat("  Warning: No SNPs available in outcome", out_name, "Skipping...\n")
      next
    }
    
    out_data$outcome <- out_name
    
    # Step 3: Harmonise Data
    dat <- harmonise_data(
      exposure_dat = exp_data, 
      outcome_dat = out_data
    )
    
    if (nrow(dat) == 0) {
      cat("  Warning: No SNPs remained after harmonisation.\n")
      next
    }
    
    # Step 4: Perform Mendelian Randomization
    res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
    
    if (nrow(res) > 0) {
      # Calculate OR and CI
      res$OR <- exp(res$b)
      res$OR_lci <- exp(res$b - 1.96 * res$se)
      res$OR_uci <- exp(res$b + 1.96 * res$se)
      
      # Step 5: Heterogeneity & Pleiotropy
      het <- mr_heterogeneity(dat)
      pleio <- mr_pleiotropy_test(dat)
      
      # Bind to internal results table for this pair
      res$pleiotropy_pval <- ifelse(nrow(pleio)>0, pleio$pval, NA)
      res$heterogeneity_pval <- ifelse(nrow(het[het$method == "Inverse variance weighted",])>0, 
                                       het[het$method == "Inverse variance weighted",]$Q_pval, NA)
      
      all_results[[paste(exp_name, out_name, sep="_")]] <- res
      
      cat("  --> IVW OR:", exp(res[res$method=="Inverse variance weighted", "b"]), 
          "| P-val:", res[res$method=="Inverse variance weighted", "pval"], "\n")
    }
  }
}

# Compile and Save
if (length(all_results) > 0) {
  final_df <- rbindlist(all_results, fill = TRUE)
  out_path <- file.path(out_folder, "GWAS_Mendelian_Randomization_Results.csv")
  fwrite(final_df, out_path)
  cat("\n=======================================================\n")
  cat("MR Pipeline Complete! Results saved to", out_path, "\n")
} else {
  cat("\nMR Pipeline completed, but no valid tests resulted in output.\n")
}

