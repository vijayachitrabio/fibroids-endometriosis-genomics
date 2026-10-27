# Check that rg with effective N is unchanged when heritability is converted to the liability scale
# (sample.prev = 0.5, as recommended for effective N; population prevalences are illustrative, rg does not depend on them)
suppressPackageStartupMessages(library(GenomicSEM))
R  <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids"
LD <- file.path(R, "gwas_scripts/eur_w_ld_chr/")
IN <- file.path(R, "Funciton_int_ge/revision_reruns/ldsc")
setwd(file.path(IN, "genomicsem_corrected_N"))
out <- list()
for (pp in list(c(0.10, 0.10), c(0.05, 0.25))) {
  fit <- ldsc(traits = file.path(IN, c("ENDO_neff.sumstats.gz", "UF_neff.sumstats.gz")),
              sample.prev = c(0.5, 0.5), population.prev = pp, ld = LD, wld = LD,
              trait.names = c("Endometriosis", "Uterine_Fibroids"), ldsc.log = paste0("genomicsem_neff_liability_", paste(pp, collapse = "_")))
  S <- fit$S; V <- fit$V; ratio <- tcrossprod(1 / sqrt(diag(S)))
  out[[length(out) + 1]] <- data.frame(pop_prev_endo = pp[1], pop_prev_uf = pp[2], rg = (S * ratio)[2, 1],
                                       rg_se = sqrt(V[2, 2]) * ratio[2, 1], h2_liab_endo = S[1, 1], h2_liab_uf = S[2, 2])
}
res <- do.call(rbind, out); print(res, digits = 6)
write.csv(res, "GenomicSEM_neff_liability_scale_check.csv", row.names = FALSE)
