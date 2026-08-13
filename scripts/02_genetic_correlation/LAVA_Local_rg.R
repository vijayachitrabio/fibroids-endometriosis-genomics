###############################################################################
# LAVA Local Genetic Correlation Analysis
# Replaces: "exploratory locus-window cross-trait correlation analysis"
#
# PURPOSE  : Compute formal local genetic correlations between Endometriosis
#            and Uterine Fibroids at 137 lead GWAS loci using the LAVA package.
#            LAVA uses LDSC-derived bivariate regression within user-defined
#            genomic windows, producing local rg with SEs and P-values that
#            are directly comparable to published literature.
#
# METHOD   : LAVA (Local Analysis of [co]Variant Associations)
#            Werme et al. (2022) Nature Genetics doi:10.1038/s41588-022-01017-y
#
# INPUTS   : FinnGen R9 GWAS: N14_ENDOMETRIOSIS.gz, CD2_BENIGN_LEIOMYOMA_UTERI.gz
#            Lead loci from coloc pipeline (outputs_v2/Regional_Correlation/
#                                           Regional_Genetic_Correlation_Results.csv)
#            LD reference: gwas_scripts/eur_w_ld_chr/ (HapMap3 EUR LD scores)
#
# OUTPUTS  : outputs_v2/LAVA_Local_rg_Results.csv  — per-locus rg, SE, P
#            outputs_v2/Regional_Correlation/Figure_LAVA_Manhattan.png
#            outputs_v2/Regional_Correlation/Figure_LAVA_Top20_loci.png
#
# REQUIRES : LAVA, data.table, ggplot2, dplyr
#
# INSTALL  : devtools::install_github("josefin-werme/LAVA")
# USAGE    : Rscript LAVA_Local_rg.R
###############################################################################

# ── Install LAVA if not present ───────────────────────────────────────────────
if (!requireNamespace("LAVA", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("josefin-werme/LAVA")
}

library(LAVA)
library(data.table)
library(ggplot2)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  BASE_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    BASE_DIR <- dirname(normalizePath(sub("--file=", "", file_arg[1])))
  } else {
    BASE_DIR <- getwd()
  }
}
GWAS_DIR  <- BASE_DIR
OUT_DIR   <- file.path(dirname(BASE_DIR), "outputs_v2")
REG_DIR   <- file.path(OUT_DIR, "Regional_Correlation")
LD_DIR    <- file.path(GWAS_DIR, "eur_w_ld_chr")
dir.create(REG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Sample sizes (FinnGen R9) ─────────────────────────────────────────────────
N_ENDO  <- 239335
N_FIBS  <- 218728

# ── Step 1: Prepare sumstats in LAVA format ───────────────────────────────────
# LAVA expects: SNP, CHR, BP, A1, A2, N, Z (or BETA + SE)
prepare_lava_sumstats <- function(gwas_file, n_total, phenotype_name) {
  cat("Preparing", phenotype_name, "...\n")
  dt <- fread(gwas_file, sep = "\t", header = TRUE)
  setnames(dt,
    old = c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval"),
    new = c("CHR","BP","A2","A1","SNP","BETA","SE","P"),
    skip_absent = TRUE)
  dt[, SNP := tstrsplit(SNP, ",", fixed = TRUE)[[1]]]
  dt[, CHR  := as.integer(CHR)]
  dt[, BP   := as.integer(BP)]
  dt[, BETA := as.numeric(BETA)]
  dt[, SE   := as.numeric(SE)]
  dt[, Z    := BETA / SE]
  dt[, N    := n_total]
  dt <- dt[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(SNP) & SNP != ""]
  out <- dt[, .(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P)]

  # Write to temp file for LAVA to read
  tmp_path <- file.path(tempdir(), paste0(phenotype_name, "_lava.txt.gz"))
  fwrite(out, tmp_path, sep = "\t", compress = "gzip")
  cat("  →", nrow(out), "SNPs written to", tmp_path, "\n")
  return(tmp_path)
}

endo_file <- file.path(GWAS_DIR, "finngen_R9_N14_ENDOMETRIOSIS.gz")
fibs_file <- file.path(GWAS_DIR, "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz")

endo_path <- prepare_lava_sumstats(endo_file, N_ENDO, "Endometriosis")
fibs_path <- prepare_lava_sumstats(fibs_file, N_FIBS, "Uterine_Fibroids")

# ── Step 2: Define locus windows ─────────────────────────────────────────────
# Load the 137 lead loci from the coloc/regional correlation analysis
loci_source <- file.path(OUT_DIR, "Regional_Correlation", "Regional_Genetic_Correlation_Results.csv")
if (file.exists(loci_source)) {
  lead_loci <- read.csv(loci_source)
  cat("Loaded", nrow(lead_loci), "lead loci from existing analysis.\n")
} else {
  stop("Lead loci file not found: ", loci_source)
}

# Build locus input for LAVA: LOC, CHR, START, STOP
WINDOW_KB <- 500    # ±500 kb around each lead SNP
locus_df <- data.frame(
  LOC   = lead_loci$locus_id,
  CHR   = lead_loci$chrom,
  START = pmax(1, lead_loci$pos - WINDOW_KB * 1000),
  STOP  = lead_loci$pos + WINDOW_KB * 1000
)
locus_file <- file.path(tempdir(), "lava_loci.txt")
write.table(locus_df, locus_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Defined", nrow(locus_df), "locus windows (±", WINDOW_KB, "kb)\n")

# ── Step 3: Load LAVA input object ───────────────────────────────────────────
input <- process.input(
  input.info.file = NULL,       # provide inline
  sample.overlap.file = NULL,   # no sample overlap between FinnGen phenotypes
  ref.prefix = file.path(LD_DIR, ""),
  phenos = c("Endometriosis", "Uterine_Fibroids")
)

# ── Step 4: Run LAVA locus-by-locus ──────────────────────────────────────────
loci <- read.loci(locus_file)
n_loc <- nrow(loci)
cat("Running LAVA on", n_loc, "loci...\n")

# Bonferroni threshold
p_bonf <- 0.05 / n_loc
cat("Bonferroni threshold: P <", signif(p_bonf, 3), "\n")

results_lava <- vector("list", n_loc)

for (i in seq_len(n_loc)) {
  loc <- loci[i, ]
  tryCatch({
    # Univariate test first (filter loci with insufficient local heritability)
    univ <- run.univ(locus = loc, input = input)
    # Only run bivariate if both phenotypes have h2 > 0 at this locus
    if (all(univ$h2 > 0, na.rm = TRUE)) {
      biv <- run.bivar(locus = loc, input = input)
      results_lava[[i]] <- data.frame(
        locus_id     = as.character(loc$LOC),
        chrom        = loc$CHR,
        start        = loc$START,
        stop         = loc$STOP,
        n_snps       = biv$n.snps,
        local_rg     = round(biv$rho, 4),
        local_rg_SE  = round(biv$rho.se, 4),
        p_val        = biv$p,
        h2_endo      = round(univ$h2[1], 6),
        h2_fibs      = round(univ$h2[2], 6),
        significant_bonferroni = biv$p < p_bonf,
        stringsAsFactors = FALSE
      )
      if (i %% 20 == 0) cat("  Processed", i, "/", n_loc, "loci\n")
    }
  }, error = function(e) {
    cat("  ✗ Locus", i, "error:", conditionMessage(e), "\n")
  })
}

# ── Step 5: Collate and annotate ──────────────────────────────────────────────
results_df <- do.call(rbind, Filter(Negate(is.null), results_lava))
results_df <- results_df[order(results_df$p_val), ]

# Add gene annotation from lead loci table
if ("locus_rsid" %in% names(lead_loci)) {
  results_df <- merge(results_df, lead_loci[, c("locus_id","locus_rsid")],
                      by = "locus_id", all.x = TRUE)
}

# FDR correction
results_df$fdr_p <- p.adjust(results_df$p_val, method = "BH")

write.csv(results_df, file.path(OUT_DIR, "LAVA_Local_rg_Results.csv"), row.names = FALSE)
cat("\n✓ LAVA results saved to:", file.path(OUT_DIR, "LAVA_Local_rg_Results.csv"), "\n")
cat("  Significant loci (Bonferroni):", sum(results_df$significant_bonferroni, na.rm = TRUE), "\n")
cat("  Significant loci (FDR < 0.05):", sum(results_df$fdr_p < 0.05, na.rm = TRUE), "\n")

# ── Step 6: Plots ─────────────────────────────────────────────────────────────
# Manhattan-style plot of local rg across the genome
results_df$chrom_num <- as.integer(results_df$chrom)
results_df <- results_df %>%
  mutate(chrom_f = factor(chrom_num)) %>%
  arrange(chrom_num, start)

# Cumulative BP positions for x-axis
chrom_offsets <- results_df %>%
  group_by(chrom_num) %>%
  summarise(max_pos = max(stop)) %>%
  mutate(offset = cumsum(lag(max_pos, default = 0))) %>%
  select(chrom_num, offset)
results_df <- merge(results_df, chrom_offsets, by = "chrom_num")
results_df$plot_pos <- results_df$start + results_df$offset

p_manhattan <- ggplot(results_df, aes(x = plot_pos, y = local_rg, colour = chrom_f)) +
  geom_point(aes(size = -log10(p_val + 1e-300)), alpha = 0.8) +
  geom_errorbar(aes(ymin = local_rg - 1.96 * local_rg_SE,
                    ymax = local_rg + 1.96 * local_rg_SE), width = 0, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", colour = "grey70") +
  scale_colour_manual(values = rep(c("#1565C0","#C62828"), 12), guide = "none") +
  scale_size_continuous(range = c(1.5, 5), name = "-log10(P)") +
  labs(title = "LAVA Local Genetic Correlation: Endometriosis vs Uterine Fibroids",
       subtitle = paste0("137 lead loci, ±500 kb windows. Bonferroni threshold P<",
                         signif(p_bonf, 2)),
       x = "Genomic position (by chromosome)", y = expression(Local~r[g])) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(REG_DIR, "Figure_LAVA_Manhattan.png"),
       p_manhattan, width = 14, height = 6, dpi = 150)

# Top 20 loci bar chart
top20 <- results_df %>% filter(!is.na(local_rg)) %>%
  top_n(-20, p_val) %>% arrange(local_rg)

p_top20 <- ggplot(top20, aes(x = reorder(locus_id, local_rg), y = local_rg,
                              fill = local_rg > 0)) +
  geom_col() +
  geom_errorbar(aes(ymin = local_rg - 1.96 * local_rg_SE,
                    ymax = local_rg + 1.96 * local_rg_SE), width = 0.4) +
  scale_fill_manual(values = c("#C62828","#1565C0"),
                    labels = c("Negative", "Positive"), name = "Direction") +
  coord_flip() +
  labs(title = "Top 20 loci by LAVA local genetic correlation P-value",
       subtitle = "Endometriosis vs Uterine Fibroids",
       x = NULL, y = expression(Local~r[g])) +
  theme_bw(base_size = 12)
ggsave(file.path(REG_DIR, "Figure_LAVA_Top20_loci.png"),
       p_top20, width = 10, height = 8, dpi = 150)

cat("✓ LAVA plots saved to:", REG_DIR, "\n")
cat("\nTop 10 loci by significance:\n")
print(head(results_df[, c("locus_id","local_rg","local_rg_SE","p_val","fdr_p",
                            "significant_bonferroni")], 10))
