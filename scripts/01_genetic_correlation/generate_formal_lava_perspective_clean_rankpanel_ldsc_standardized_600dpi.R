suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

base_dir <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/uterine_fibroids"
lava_dir <- file.path(base_dir, "outputs_v2", "formal_ld_audit", "lava")
pub_dir <- file.path(base_dir, "outputs_v2", "Figures_Publication", "Figure3_LAVA_Panels")
asset_dir <- file.path(base_dir, "manuscript_assets_latest_2026_05_20", "figures_main")

dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

lava <- read.csv(
  file.path(lava_dir, "Table_LAVA_Full_2495Blocks_AllResults_overlap_ldsc_standardized_withFDR.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
loci <- read.delim(
  file.path(lava_dir, "lava_loci_2495_blocks.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

for (nm in c("CHR", "START", "STOP", "n_snps", "K", "rho", "rho.lower",
             "rho.upper", "r2", "r2.lower", "r2.upper", "p", "p_fdr")) {
  if (nm %in% names(lava)) lava[[nm]] <- as.numeric(lava[[nm]])
}
for (nm in c("CHR", "START", "STOP")) loci[[nm]] <- as.numeric(loci[[nm]])

chr_lengths <- aggregate(STOP ~ CHR, loci, max)
chr_lengths <- chr_lengths[order(chr_lengths$CHR), ]
chr_lengths$offset <- c(0, head(cumsum(chr_lengths$STOP), -1))
chr_lengths$center <- chr_lengths$offset + chr_lengths$STOP / 2

lava$mid <- (lava$START + lava$STOP) / 2
lava <- merge(lava, chr_lengths[, c("CHR", "offset")], by = "CHR", all.x = TRUE)
lava$x <- lava$mid + lava$offset
lava$sig <- !is.na(lava$p_fdr) & lava$p_fdr < 0.05
lava$rho.lower.plot <- pmax(lava$rho.lower, -1)
lava$rho.upper.plot <- pmin(lava$rho.upper, 1)
lava$neglog10p <- -log10(lava$p)

lava$locus <- ""
lava$locus[lava$CHR == 1 & lava$START <= 22750000 & lava$STOP >= 22100000] <- "WNT4"
lava$locus[lava$CHR == 6 & lava$START <= 153800000 & lava$STOP >= 151600000] <- "ESR1"
lava$locus[lava$CHR == 11 & lava$START <= 30600000 & lava$STOP >= 29500000] <- "WT1"

lava$label <- ""
for (lab in unique(lava$locus[lava$locus != ""])) {
  idx <- which(lava$locus == lab)
  lava$label[idx[which.min(lava$p[idx])]] <- lab
}

bonf_line <- -log10(0.05 / nrow(loci))
sig_col <- "#2166ac"

theme_clean <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", color = "black", size = base_size + 1),
      plot.title.position = "plot",
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
}

panel_a <- ggplot(lava, aes(x = x, y = rho)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_linerange(aes(ymin = rho.lower.plot, ymax = rho.upper.plot),
                 color = "grey78", linewidth = 0.22, alpha = 0.85) +
  geom_point(aes(size = sig, fill = sig), shape = 21, color = "black",
             stroke = 0.28, alpha = 0.95) +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 3.0)) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = sig_col)) +
  scale_x_continuous(
    breaks = chr_lengths$center,
    labels = chr_lengths$CHR,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(limits = c(-1.05, 1.05), breaks = seq(-1, 1, 0.5)) +
  labs(title = "A  Local rg by chromosome",
       x = "Chromosome",
       y = "Local genetic correlation (rho)") +
  theme_clean(11)

lava_ranked <- lava[order(lava$rho, lava$p), ]
lava_ranked$rank <- seq_len(nrow(lava_ranked))

panel_b <- ggplot(lava_ranked, aes(x = rank, y = rho)) +
  geom_hline(yintercept = median(lava$rho, na.rm = TRUE),
             color = "grey35", linetype = "dashed", linewidth = 0.35) +
  geom_segment(aes(xend = rank, y = 0, yend = rho),
               color = "grey82", linewidth = 0.25) +
  geom_point(
    aes(size = sig, fill = sig),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  scale_size_manual(values = c("FALSE" = 2.1, "TRUE" = 3.1)) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = sig_col)) +
  scale_x_continuous(breaks = c(1, 10, 20, nrow(lava_ranked))) +
  scale_y_continuous(limits = c(-1.05, 1.05), breaks = seq(-1, 1, 0.5)) +
  labs(title = "B  Ranked local rg",
       x = "LAVA blocks ranked by rho",
       y = "Local genetic correlation (rho)") +
  theme_clean(11)

panel_c <- ggplot(lava, aes(x = rho, y = neglog10p)) +
  geom_hline(yintercept = bonf_line, color = "grey45", linetype = "dashed", linewidth = 0.35) +
  geom_point(aes(size = sig, fill = sig), shape = 21, color = "black",
             stroke = 0.28, alpha = 0.95) +
  geom_text_repel(
    data = subset(lava, label != ""),
    aes(label = label),
    color = "black",
    size = 3.3,
    min.segment.length = 0,
    segment.color = "grey45",
    segment.size = 0.25,
    box.padding = 0.25,
    max.overlaps = Inf,
    seed = 42
  ) +
  scale_size_manual(values = c("FALSE" = 2.0, "TRUE" = 3.2)) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = sig_col)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(xlim = c(-1.02, 1.02)) +
  labs(title = "C  LAVA significance",
       x = "Local genetic correlation (rho)",
       y = expression(-log[10](italic(P)))) +
  theme_clean(11)

combined <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(1.08, 1))

save_plot <- function(plot, stem, width, height) {
  for (out_dir in c(pub_dir, asset_dir)) {
    ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
           width = width, height = height, units = "in",
           dpi = 600, bg = "white", limitsize = FALSE)
    ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
           width = width, height = height, units = "in",
           bg = "white", limitsize = FALSE)
  }
}

save_plot(combined, "Figure_3_LAVA_Local_rg_LDSC_Standardized_RANKPANEL_600dpi", 15.3, 9.05)
save_plot(panel_a, "Figure_3A_LAVA_Local_rg_Chromosomal_LDSC_Standardized_600dpi", 10.5, 4.56)
save_plot(panel_b, "Figure_3B_LAVA_Local_rg_Ranked_LDSC_Standardized_600dpi", 5.95, 4.56)
save_plot(panel_c, "Figure_3C_LAVA_Local_rg_Significance_LDSC_Standardized_600dpi", 6.5, 4.87)

message("Saved LDSC-standardized overlap-aware LAVA rank-panel 600 dpi figure set.")
