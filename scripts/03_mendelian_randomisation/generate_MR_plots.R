library(TwoSampleMR)
library(data.table)
library(ggplot2)

# OPENGWAS_JWT loaded from ~/.Renviron — see MR_phenotype_divergence.R for setup
Sys.setenv(OPENGWAS_JWT = Sys.getenv("OPENGWAS_JWT"))
if (nchar(Sys.getenv("OPENGWAS_JWT")) == 0)
  stop("OPENGWAS_JWT is not set. Add it to ~/.Renviron and restart R.")

cat("Starting MR Plot Generation Pipeline...\n")

# Same file names as our local pipeline
exp_files <- list(
  Endometriosis = "finngen_R9_N14_ENDOMETRIOSIS.gz",
  Uterine_Fibroids = "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz"
)

out_files <- list(
  IBS = "finngen_R9_K11_IBS.gz",
  Endometrial_Cancer = "finngen_R9_C3_CORPUS_UTERI_EXALLC.gz"
)

out_dir <- "../outputs_v2/MR_Plots"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Function to safely generate and save plots
generate_all_mr_plots <- function(dat, res, exp_name, out_name) {
  cat("  Generating plots for", exp_name, "->", out_name, "\n")
  
  prefix <- file.path(out_dir, paste0(exp_name, "_vs_", out_name))
  
  # 1. Scatter Plot
  p1 <- mr_scatter_plot(res, dat)
  ggsave(paste0(prefix, "_scatter.png"), plot = p1[[1]], width = 7, height = 7)
  
  # 2. Single SNP Forest Plot
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)
  ggsave(paste0(prefix, "_forest.png"), plot = p2[[1]], width = 7, height = 9)
  
  # 3. Funnel Plot
  p3 <- mr_funnel_plot(res_single)
  ggsave(paste0(prefix, "_funnel.png"), plot = p3[[1]], width = 7, height = 7)
  
  # 4. Leave-one-out Plot
  res_loo <- mr_leaveoneout(dat)
  p4 <- mr_leaveoneout_plot(res_loo)
  ggsave(paste0(prefix, "_leave_one_out.png"), plot = p4[[1]], width = 7, height = 9)
}

# Run the rapid extraction loop to get the 'dat' objects again
for (exp_name in names(exp_files)) {
  exp_file <- exp_files[[exp_name]]
  
  cat("\nLoading", exp_name, "...\n")
  d_exp <- fread(exp_file)
  sig_exp <- d_exp[pval < 5e-8]
  
  exp_dat <- format_data(
    as.data.frame(sig_exp),
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "phenotype",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    samplesize_col = "N",
    eaf_col = "af_alt"
  )
  exp_dat$id.exposure <- exp_name
  exp_dat$exposure <- exp_name
  
  cat("  Clumping locally via API...\n")
  clumped_exp <- clump_data(
    exp_dat,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR"
  )
  
  if (nrow(clumped_exp) == 0) next
  clumped_snps <- clumped_exp$SNP
  
  for (out_name in names(out_files)) {
    out_file <- out_files[[out_name]]
    cat("  Against", out_name, "...\n")
    
    d_out <- fread(out_file)
    sig_out <- d_out[rsids %in% clumped_snps]
    
    if (nrow(sig_out) == 0) next
    
    out_dat <- format_data(
      as.data.frame(sig_out),
      type = "outcome",
      snps = clumped_snps,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "sebeta",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval",
      eaf_col = "af_alt"
    )
    out_dat$id.outcome <- out_name
    out_dat$outcome <- out_name
    
    dat <- harmonise_data(exposure_dat = clumped_exp, outcome_dat = out_dat)
    
    # Run primary MR result
    res <- mr(dat)
    
    # Send dat and res to plotting function
    generate_all_mr_plots(dat, res, exp_name, out_name)
  }
}

cat("\nPlotting Complete! All images saved to outputs_v2/MR_Plots/\n")
