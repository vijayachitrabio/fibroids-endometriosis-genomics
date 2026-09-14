library(data.table)
library(GenomicSEM)

cat("===================================================\n")
cat("Starting Global Genetic Correlation (LDSC) Pipeline\n")
cat("===================================================\n")

out_dir <- "../outputs_v2"
ld_folder <- "eur_w_ld_chr"

# 1. Clean FinnGen Summary Statistics for GenomicSEM Munging
phenos <- list(
  Endometriosis = "finngen_R9_N14_ENDOMETRIOSIS.gz",
  Uterine_Fibroids = "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz"
)

# FinnGen Total Sample Sizes (from R9 Manifest)
n_totals <- list(
  Endometriosis = 239335,
  Uterine_Fibroids = 218728
)

cleaned_files <- c()
for (pName in names(phenos)) {
  file_in <- phenos[[pName]]
  file_out <- file.path(out_dir, paste0("Cleaned_", pName, ".txt"))
  
  if (!file.exists(file_out)) {
    cat(sprintf("Cleaning raw 3GB arrays for %s...\n", pName))
    dt <- fread(file_in, select = c("rsids", "alt", "ref", "beta", "sebeta", "pval", "af_alt"))
    
    # Standardize column names for munge mapping
    setnames(dt,
             old = c("rsids", "alt", "ref", "beta", "sebeta", "pval", "af_alt"),
             new = c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF"))
    
    # Strip empty rsids
    dt <- dt[SNP != "" & !is.na(SNP) & !is.na(P)]
    
    # Add absolute Sample Size column required for Heritability scaling
    dt[, N := n_totals[[pName]]]
    
    fwrite(dt, file_out, sep="\t")
  }
  cleaned_files <- c(cleaned_files, file_out)
}

cat("\nExecuting native LDSC Munging against 1000G Linkage Maps...\n")
sumstats_files <- c()
for (i in seq_along(names(phenos))) {
  pName <- names(phenos)[i]
  m_out <- file.path(out_dir, paste0(pName, "_munged"))
  sumstats_files <- c(sumstats_files, paste0(m_out, ".sumstats.gz"))
  
  if (!file.exists(paste0(m_out, ".sumstats.gz"))) {
    munge(files = cleaned_files[i], 
          hm3 = file.path(ld_folder, "w_hm3.snplist"), 
          trait.names = pName, 
          N = n_totals[[pName]], 
          info.filter = 0.9, 
          maf.filter = 0.01)
          
    # Munge naturally dumps in current directory, so we move it to outputs manually
    file.rename(paste0(pName, ".sumstats.gz"), paste0(m_out, ".sumstats.gz"))
  }
}

cat("\nCalculating Global Genetic Correlation (rg)...\n")

# Run actual cross-trait LDSC 
ldsc_res <- tryCatch({
  ldsc(traits = sumstats_files, 
       sample.prev = c(NA, NA), 
       population.prev = c(NA, NA), 
       ld = file.path(ld_folder), 
       wld = file.path(ld_folder), 
       trait.names = names(phenos))
}, error = function(e){
  cat("LDSC Execution Error:", e$message, "\n")
  NULL
})

if (!is.null(ldsc_res)) {
  cat("\n===================================================\n")
  cat("LDSC Output Successful!\n")
  print(ldsc_res$S)
  print(ldsc_res$V)
  
  # NOTE: ldsc_res$S is the genetic COVARIANCE matrix, not the correlation matrix.
  # rg = S[1,2] / sqrt(S[1,1] * S[2,2])  (calculated below)
  write.csv(ldsc_res$S, file.path(out_dir, "LDSC_Genetic_Covariance_Matrix.csv"), row.names = TRUE)
  write.csv(ldsc_res$V, file.path(out_dir, "LDSC_Correlation_Variance_SE.csv"), row.names = TRUE)
  
  # The genetic correlation (rg) is usually S[1,2] / sqrt(S[1,1] * S[2,2])
  # But GenomicSEM provides it inherently if we check its structural objects
  rg <- as.numeric((ldsc_res$S[1,2]) / sqrt(ldsc_res$S[1,1] * ldsc_res$S[2,2]))
  cat("\nCalculated Global Genetic Correlation (rg):", round(rg, 4), "\n")
  
} else {
  cat("\nFailure generating LDSC architectures.\n")
}
