# GenomicSEM LDSC with the munged-file, total (R9) and effective (R9) sample sizes
suppressPackageStartupMessages(library(GenomicSEM))
R  <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids"
LD <- file.path(R, "gwas_scripts/eur_w_ld_chr/")
IN <- file.path(R, "Funciton_int_ge/revision_reruns/ldsc")
OUT <- file.path(IN, "genomicsem_corrected_N")
setwd(OUT)
res <- list()
for (s in c("submitted", "total", "neff")) {
  fit <- ldsc(traits = file.path(IN, c(paste0("ENDO_", s, ".sumstats.gz"), paste0("UF_", s, ".sumstats.gz"))),
              sample.prev = c(NA, NA), population.prev = c(NA, NA), ld = LD, wld = LD,
              trait.names = c("Endometriosis", "Uterine_Fibroids"), ldsc.log = paste0("genomicsem_", s))
  S <- fit$S; V <- fit$V
  ratio <- tcrossprod(1 / sqrt(diag(S)))
  rg <- (S * ratio)[2, 1]
  se_rg <- sqrt(V[2, 2]) * ratio[2, 1]
  se_cov <- sqrt(V[2, 2])
  res[[s]] <- data.frame(sample_size = s, gcov = S[2, 1], gcov_se = se_cov, gcov_p = 2 * pnorm(-abs(S[2, 1] / se_cov)),
                         rg = rg, rg_se = se_rg, h2_endo = S[1, 1], h2_uf = S[2, 2],
                         intercept_endo = fit$I[1, 1], intercept_uf = fit$I[2, 2], intercept_cross = fit$I[2, 1])
}
out <- do.call(rbind, res)
print(out, digits = 4)
write.csv(out, file.path(OUT, "GenomicSEM_rg_by_sample_size.csv"), row.names = FALSE)
