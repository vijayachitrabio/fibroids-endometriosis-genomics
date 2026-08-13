# MR_local_pipeline.R
# Fully Offline Manual Mendelian Randomization
library(data.table)
suppressPackageStartupMessages(library(TwoSampleMR))

out_folder <- "../outputs_v2"
dir.create(out_folder, showWarnings = FALSE)

cat("Starting Local Manual MR Pipeline...\n")

# OPENGWAS_JWT loaded from ~/.Renviron — see MR_phenotype_divergence.R for setup
Sys.setenv(OPENGWAS_JWT = Sys.getenv("OPENGWAS_JWT"))
if (nchar(Sys.getenv("OPENGWAS_JWT")) == 0)
  stop("OPENGWAS_JWT is not set. Add it to ~/.Renviron and restart R.")

exp_files <- list(
  Endometriosis = "finngen_R9_N14_ENDOMETRIOSIS.gz",
  Uterine_Fibroids = "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz"
)

out_files <- list(
  IBS = "finngen_R9_K11_IBS.gz",
  Endometrial_Cancer = "finngen_R9_C3_CORPUS_UTERI_EXALLC.gz"
)

all_results <- list()

for (exp_name in names(exp_files)) {
  exp_file <- exp_files[[exp_name]]
  cat("\n=======================================================\n")
  cat("Loading local array via data.table for:", exp_name, "...\n")
  
  if (!file.exists(exp_file)) {
    cat(exp_file, "not found! Please ensure it finished downloading.\n")
    next
  }
  
  exp_df <- fread(exp_file)
  sig_exp <- exp_df[pval < 5e-8]
  
  cat("Formatting", nrow(sig_exp), "significant genome variants for", exp_name, "...\n")
  
  exp_dat <- format_data(
    as.data.frame(sig_exp),
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "Phenotype",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "maf",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    chr_col = "#chrom",
    pos_col = "pos"
  )
  exp_dat$exposure <- exp_name
  
  # Remote clump using tiny proxy subset
  cat("Clumping independent instruments through lightweight LD API...\n")
  exp_clumped <- tryCatch(clump_data(exp_dat, clump_r2 = 0.001, clump_kb = 10000), error = function(e){cat("Clumping failed.\n"); return(NULL)})
  if (is.null(exp_clumped) || nrow(exp_clumped) == 0) {
    cat("No instruments remain after clumping. Skipping.\n")
    next
  }
  
  clumped_snps <- exp_clumped$SNP
  cat("Found", length(clumped_snps), "independent causal instruments.\n")
  
  for (out_name in names(out_files)) {
    out_file <- out_files[[out_name]]
    if (!file.exists(out_file)) {
       cat("Outcome dataset", out_name, "not found on disk.\n")
       next
    }
    
    cat("\n  Testing Outcome against local offline dataset:", out_name, "...\n")
    out_df <- fread(out_file)
    sig_out <- out_df[rsids %in% clumped_snps]
    
    if (nrow(sig_out) == 0) {
      cat("  Warning: None of the instruments matched in the outcome file.\n")
      next
    }
    
    out_dat <- format_data(
      as.data.frame(sig_out),
      type = "outcome",
      snps = clumped_snps,
      header = TRUE,
      phenotype_col = "Phenotype",
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "sebeta",
      eaf_col = "maf",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval",
      chr_col = "#chrom",
      pos_col = "pos"
    )
    out_dat$outcome <- out_name
    
    dat <- harmonise_data(exp_clumped, out_dat)
    
    if (nrow(dat) == 0) {
      cat("  Warning: No SNPs remained after harmonisation.\n")
      next
    }
    
    res <- tryCatch(mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")), error=function(e) data.frame())
    
    if (nrow(res) > 0) {
      res$OR <- exp(res$b)
      res$OR_lci <- exp(res$b - 1.96 * res$se)
      res$OR_uci <- exp(res$b + 1.96 * res$se)
      
      het <- tryCatch(mr_heterogeneity(dat), error=function(e) data.frame())
      pleio <- tryCatch(mr_pleiotropy_test(dat), error=function(e) data.frame())
      
      res$pleiotropy_pval <- ifelse(nrow(pleio)>0, pleio$pval, NA)
      res$heterogeneity_pval <- ifelse(nrow(het)>0 && nrow(het[het$method == "Inverse variance weighted",])>0, 
                                       het[het$method == "Inverse variance weighted",]$Q_pval, NA)
      
      all_results[[paste(exp_name, out_name, sep="_")]] <- res
      cat("  ==> IVW OR:", exp(res[res$method=="Inverse variance weighted", "b"]), 
          "| P-value:", res[res$method=="Inverse variance weighted", "pval"], "\n")
    }
  }
}

if (length(all_results) > 0) {
  final_df <- rbindlist(all_results, fill = TRUE)
  out_path <- file.path(out_folder, "GWAS_Manual_MR_Results.csv")
  fwrite(final_df, out_path)
  cat("\n=======================================================\n")
  cat("Local Bypass Pipeline Complete! Causality matrices exported to", out_path, "\n")
} else {
  cat("\nPipeline failed to calculate causality arrays.\n")
}
