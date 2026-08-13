#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(patchwork)
})

root <- normalizePath(".", mustWork = TRUE)
out_dir <- file.path(root, "figures")
pub_dir <- file.path(root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

endo <- read_csv(file.path(root, "results", "Olink_Endo_Nominal_p001.csv"), show_col_types = FALSE) |>
  transmute(
    protein = toupper(protein),
    logfc = logfc_endo,
    p = p_endo,
    fdr = fdr_endo,
    panel = "Endometriosis"
  )

fib <- read_csv(file.path(root, "results", "Olink_Fibroids_Nominal_p001.csv"), show_col_types = FALSE) |>
  transmute(
    protein = toupper(protein),
    logfc = logfc_fib,
    p = p_fib,
    fdr = fdr_fib,
    panel = "Uterine fibroids"
  )

plot_data <- bind_rows(endo, fib) |>
  mutate(
    panel = factor(panel, levels = c("Endometriosis", "Uterine fibroids")),
    neg_log10_p = -log10(p),
    status = case_when(
      fdr < 0.05 & logfc > 0 ~ "Increased",
      fdr < 0.05 & logfc < 0 ~ "Decreased",
      TRUE ~ "Nominal"
    ),
    status = factor(status, levels = c("Increased", "Decreased", "Nominal")),
    label = if_else(fdr < 0.05, protein, NA_character_)
  )

p_threshold <- 0.001

p <- ggplot(plot_data, aes(x = logfc, y = neg_log10_p)) +
  geom_hline(yintercept = -log10(p_threshold), linewidth = 0.35, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey65") +
  geom_point(
    aes(fill = status, size = status),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.95
  ) +
  geom_text_repel(
    aes(label = label),
    color = "black",
    size = 2.15,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.25,
    box.padding = 0.32,
    point.padding = 0.18,
    force = 12,
    max.iter = 10000,
    max.time = 2,
    max.overlaps = Inf,
    seed = 42,
    na.rm = TRUE
  ) +
  facet_wrap(~panel, nrow = 1, scales = "free") +
  scale_fill_manual(
    values = c("Increased" = "#D73027", "Decreased" = "#4575B4", "Nominal" = "grey85"),
    breaks = c("Increased", "Decreased", "Nominal")
  ) +
  scale_size_manual(
    values = c("Increased" = 2.4, "Decreased" = 2.4, "Nominal" = 1.8),
    breaks = c("Increased", "Decreased", "Nominal")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "Log fold-change", y = expression(-log[10](P))) +
  theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 11, face = "plain"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 9),
    legend.key.width = unit(0.45, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.spacing.x = unit(1.1, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 6, r = 8, b = 4, l = 8)
  )

base <- "FigS20_Olink_Volcano_MINIMAL_600dpi"
png_out <- file.path(out_dir, paste0(base, ".png"))
pdf_out <- file.path(out_dir, paste0(base, ".pdf"))
pub_png <- file.path(pub_dir, paste0(base, ".png"))
pub_pdf <- file.path(pub_dir, paste0(base, ".pdf"))

ggsave(png_out, p, width = 9.0, height = 4.8, units = "in", dpi = 600, bg = "white")
ggsave(pdf_out, p, width = 9.0, height = 4.8, units = "in", device = "pdf", bg = "white")
ggsave(pub_png, p, width = 9.0, height = 4.8, units = "in", dpi = 600, bg = "white")
ggsave(pub_pdf, p, width = 9.0, height = 4.8, units = "in", device = "pdf", bg = "white")

message("Saved: ", png_out)
message("Saved: ", pdf_out)
