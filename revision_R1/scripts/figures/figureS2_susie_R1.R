suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
  library(scales)
})

root <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids"
fig <- file.path(root, "Funciton_int_ge")
src <- file.path(root, "outputs_v2/formal_ld_audit/susie_finemap/SuSiE_Shared_Loci_Variant_PIPs.csv")
out <- file.path(fig, "revision_figures")

d <- fread(src)
# Credible-set membership from the repeat run (Supplementary Table S3, Panel F); P values from the full variant set
mem <- fread(file.path(fig, "revision_reruns/susie/SuSiE_credible_set_membership_L10.csv"))[, .(locus, trait, SNP, credible_set_rerun = credible_set)]
d[, c("credible_set", "cs_label", "pip") := NULL]
d <- merge(d, mem, by = c("locus", "trait", "SNP"), all.x = TRUE)
setnames(d, "credible_set_rerun", "credible_set")
conv <- fread(file.path(fig, "revision_reruns/susie/SuSiE_L_sensitivity_rerun.csv"))[L == 10, .(locus, trait, converged)]
d <- d[is.finite(P) & P > 0 & is.finite(BP)]
d[, neglog10p := -log10(P)]
d[, cs_group := fifelse(is.na(credible_set) | trimws(credible_set) == "", "Not in credible set", as.character(credible_set))]
cs_found <- unique(d[cs_group != "Not in credible set", cs_group])
cs_levels <- c(cs_found[order(as.integer(sub("^L", "", cs_found)))], "Not in credible set")
d[, cs_group := factor(cs_group, levels = cs_levels)]
pal <- setNames(hcl.colors(max(1, length(cs_levels) - 1), "Dark 3"), cs_levels[seq_len(max(1, length(cs_levels) - 1))])
pal <- c(pal, "Not in credible set" = "#B7B7B7")

loc_order <- c("WNT4", "GREB1", "WT1", "ESR1_SYNE1", "DNM3")
trait_order <- c("Endometriosis", "Uterine_Fibroids")
plots <- list()
for (loc in loc_order) for (tr in trait_order) {
  dd <- d[locus == loc & trait == tr]
  if (!nrow(dd)) next
  locus_label <- switch(loc,
    WT1 = "chr11 (30.2 Mb region)",
    ESR1_SYNE1 = "ESR1/SYNE1",
    loc
  )
  title_locus <- if (loc == "WT1") locus_label else paste0("<i>", locus_label, "</i>")
  trait_label <- if (tr == "Uterine_Fibroids") "Uterine fibroids" else "Endometriosis"
  lead_bp <- dd[which.min(P), BP]
  chr <- unique(dd$CHR)
  plots[[length(plots) + 1]] <- ggplot(dd, aes(BP, neglog10p, colour = cs_group)) +
    geom_vline(xintercept = lead_bp, colour = "#777777", linewidth = 0.4, linetype = "dashed") +
    geom_point(data = dd[cs_group == "Not in credible set"], size = 0.5, alpha = 0.65) +
    geom_point(data = dd[cs_group != "Not in credible set"], size = 1.4, alpha = 0.95) +
    scale_colour_manual(values = pal, breaks = cs_levels, limits = cs_levels, drop = FALSE, name = "Credible set (L = 10)") +
    scale_x_continuous(labels = label_number(big.mark = ","), breaks = pretty_breaks(4), expand = expansion(mult = 0.02)) +
    labs(title = paste0(title_locus, " fine-mapping: ", trait_label,
        if (!isTRUE(conv[locus == loc & trait == tr, converged])) " (did not converge)" else ""),
      x = paste0("Chromosome ", chr, " position (GRCh38)"), y = expression(-log[10](P))) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_markdown(size = 8.5, face = "bold", margin = margin(b = 3)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.25),
      axis.text = element_text(size = 6.5, colour = "black"),
      axis.title = element_text(size = 7.5),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.key.width = unit(0.35, "cm"),
      legend.key.height = unit(0.25, "cm"),
      plot.margin = margin(5, 6, 4, 6)
    )
}
p <- wrap_plots(plots, ncol = 2) + plot_annotation(tag_levels = "A") & theme(legend.position = "none", plot.tag = element_text(size = 10, face = "bold"))

dir.create(out, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out, "Supplementary_Figure_S2_SuSiE_R1.png"), p,
  width = 11.7, height = 12.6, units = "in", dpi = 600, bg = "white")
ggsave(file.path(out, "Supplementary_Figure_S2_SuSiE_R1.pdf"), p,
  width = 11.7, height = 12.6, units = "in", bg = "white")
