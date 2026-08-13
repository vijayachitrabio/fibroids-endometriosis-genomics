# ==============================================================================
# 04_pMR_coloc_visualization.R
# Purpose : Generate all figures and supplementary tables for the
#           proteome-MR + colocalization analysis.
#

#   Figure B — Proteome-MR forest plot
#              OR (95% CI) per FDR-significant protein–outcome pair
#              Panels per outcome, proteins on y-axis
#
#   Figure C — Evidence integration heatmap
#              Proteins (rows) × outcomes (columns)
#              Cell fill = IVW OR (log scale)
#              Cell border thickness = PP.H4 coloc evidence
#              Stars = FDR significance tier
#
#   Supp Table S_pMR — Full proteome-MR results (IVW, WM, Egger)
#   Supp Table S_coloc — Full colocalization posteriors summary
#
# Author  : Vijayachitra Modhukur  |  April 2026
# Depends : ggplot2, dplyr, readr, tidyr, ggrepel, patchwork,
#           RColorBrewer, scales
# ==============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(scales)

# ── PATHS ─────────────────────────────────────────────────────────────────────
BASE_DIR  <- Sys.getenv("FIBROID_BASE_DIR", unset = getwd())
OUT_DIR   <- Sys.getenv("FIBROID_OUT_DIR",  unset = file.path(BASE_DIR, "outputs_v2"))
PLOT_DIR  <- Sys.getenv("FIBROID_PLOT_DIR", unset = file.path(BASE_DIR, "plots"))
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
coloc   <- read_csv(file.path(OUT_DIR, "coloc_results_summary.csv"),  show_col_types = FALSE)
pMR_ivw <- read_csv(file.path(OUT_DIR, "pMR_ivw_results.csv"),        show_col_types = FALSE)
pMR_all <- read_csv(file.path(OUT_DIR, "pMR_all_results.csv"),        show_col_types = FALSE)
pMR_sig <- read_csv(file.path(OUT_DIR, "pMR_FDR_significant.csv"),    show_col_types = FALSE)

# Olink direction
fibroid_logfc <- data.frame(
  protein   = c("CHRDL2","TNFRSF11B","EFEMP1","FLT3LG","EDA2R",
                "CDH3","COL9A1","FGF23","TNFRSF17","PTK7",
                "CD55","TNFSF11","LEFTY2","LPL","CHL1","TFPI"),
  logfc_fib = c(0.148,-0.059,-0.052,-0.053,-0.054,
                0.062,-0.074, 0.118,-0.058, 0.055,
               -0.033, 0.082, 0.078,-0.057,-0.036,-0.037),
  direction = c("Up","Down","Down","Down","Down",
                "Up","Down","Up","Down","Up",
                "Down","Up","Up","Down","Down","Down")
)

THEME <- theme_classic(base_size = 12) +
  theme(
    axis.text        = element_text(colour = "black"),
    strip.background = element_rect(fill = "#EEF2F7", colour = NA),
    strip.text       = element_text(face = "bold", size = 11),
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "bottom",
    legend.key.size  = unit(0.5, "cm")
  )



# ── FIGURE B: PROTEOME-MR FOREST PLOT ─────────────────────────────────────────
# Show all three estimators for FDR-significant protein-outcome pairs
forest_dat <- pMR_all %>%
  filter(paste0(protein, "_", outcome) %in%
           paste0(pMR_sig$protein, "_", pMR_sig$outcome)) %>%
  mutate(
    method   = factor(method, levels = c("IVW","WM","Egger")),
    sig_star = case_when(
      method == "IVW" & pval < 0.001 ~ "***",
      method == "IVW" & pval < 0.01  ~ "**",
      method == "IVW" & pval < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label = paste0(protein, " -> ", outcome)
  ) %>%
  left_join(select(coloc, protein, PP.H4), by = "protein") %>%
  mutate(
    coloc_annot = case_when(
      PP.H4 >= 0.80 ~ "◆ Colocalized",
      PP.H4 >= 0.50 ~ "◇ Moderate coloc",
      TRUE          ~ ""
    )
  )

# Order proteins by IVW OR within each outcome
protein_order <- forest_dat %>%
  filter(method == "IVW") %>%
  arrange(outcome, OR) %>%
  pull(label) %>% unique()

forest_dat$label <- factor(forest_dat$label, levels = rev(protein_order))

fig_B <- ggplot(forest_dat,
                aes(x = OR, y = label,
                    colour = method, shape = method)) +
  geom_vline(xintercept = 1, linetype = "solid", colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0, linewidth = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5,
             position = position_dodge(width = 0.5)) +
  geom_text(data = filter(forest_dat, method == "IVW", sig_star != ""),
            aes(label = sig_star), hjust = -0.3,
            colour = "black", size = 4.5,
            position = position_dodge(width = 0.5)) +
  facet_wrap(~outcome, scales = "free_y", ncol = 2) +
  scale_colour_manual(
    name   = "MR estimator",
    values = c("IVW" = "#2C3E50", "WM" = "#2980B9", "Egger" = "#E74C3C")
  ) +
  scale_shape_manual(
    name   = "MR estimator",
    values = c("IVW" = 16, "WM" = 17, "Egger" = 15)
  ) +
  scale_x_log10(breaks = c(0.5, 0.75, 1.0, 1.5, 2.0, 3.0),
                labels  = c("0.50","0.75","1.00","1.50","2.00","3.00")) +
  labs(
    title    = "Proteome-MR: causal effect of fibroid proteins on comorbidity risk",
    subtitle = "OR per 1-SD increase in genetically predicted protein level (cis-pQTL instruments)\nStars = IVW FDR significance: * <0.05, ** <0.01, *** <0.001",
    x        = "Odds ratio (95% CI, log scale)",
    y        = NULL
  ) +
  THEME

ggsave(file.path(PLOT_DIR, "FigB_pMR_ForestPlot.pdf"),
       fig_B, width = 12, height = max(6, nrow(pMR_sig) * 1.2), dpi = 300)
ggsave(file.path(PLOT_DIR, "FigB_pMR_ForestPlot.png"),
       fig_B, width = 12, height = max(6, nrow(pMR_sig) * 1.2), dpi = 300)
message("Saved Figure B: Proteome-MR forest plot")

# ── FIGURE C: EVIDENCE INTEGRATION HEATMAP ────────────────────────────────────
# Rows = proteins, Columns = outcomes
# Cell fill = log(IVW OR), capped at ±0.5
# Cell outline = PP.H4 tier (thick = strong coloc)
# Text = FDR star tier
heatmap_dat <- pMR_ivw %>%
  # PP.H4 is already present from the Script 03b join — avoid duplicate column
  { if (!"PP.H4" %in% names(.))
      left_join(., select(coloc, protein, PP.H4), by = "protein")
    else . } %>%
  replace_na(list(PP.H4 = 0)) %>%
  mutate(
    log_OR      = log(OR),
    log_OR_cap  = pmax(pmin(log_OR, 0.5), -0.5),
    fdr_star    = case_when(FDR < 0.001 ~ "***",
                            FDR < 0.01  ~ "**",
                            FDR < 0.05  ~ "*",
                            TRUE        ~ ""),
    coloc_bord  = case_when(PP.H4 >= 0.80 ~ "Strong",
                            PP.H4 >= 0.50 ~ "Moderate",
                            TRUE          ~ "None")
  )

# Order proteins by average absolute log(OR) across outcomes
prot_order <- heatmap_dat %>%
  group_by(protein) %>%
  summarise(mean_abs_logOR = mean(abs(log_OR), na.rm = TRUE)) %>%
  arrange(desc(mean_abs_logOR)) %>%
  pull(protein)

out_order <- c("PCOS","Ovarian cyst","Ovarian cancer",
               "Hypertension","Depression","Coeliac disease","Breast cancer")

heatmap_dat$protein <- factor(heatmap_dat$protein, levels = rev(prot_order))
heatmap_dat$outcome <- factor(heatmap_dat$outcome, levels = out_order)

fig_C <- ggplot(heatmap_dat,
                aes(x = outcome, y = protein, fill = log_OR_cap)) +
  geom_tile(aes(colour = coloc_bord,
                linewidth = coloc_bord), show.legend = TRUE) +
  geom_text(aes(label = fdr_star), size = 5, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    name     = "log(OR)\n[IVW, capped ±0.5]",
    low      = "#2471A3",
    mid      = "white",
    high     = "#C0392B",
    midpoint = 0,
    limits   = c(-0.5, 0.5)
  ) +
  scale_colour_manual(
    name   = "Coloc evidence\n(PP.H4)",
    values = c("Strong"   = "#C0392B",
               "Moderate" = "#E67E22",
               "None"     = "grey80")
  ) +
  scale_linewidth_manual(
    name   = "Coloc evidence\n(PP.H4)",
    values = c("Strong" = 2.0, "Moderate" = 1.2, "None" = 0.3)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Proteome-MR × colocalization evidence matrix",
    subtitle = "Fill = log(OR) per 1-SD protein increase | Border = pQTL-GWAS colocalization strength\nStars: * FDR<0.05, ** FDR<0.01, *** FDR<0.001",
    x        = NULL,
    y        = "Fibroid-associated protein"
  ) +
  THEME +
  theme(
    axis.text.x       = element_text(angle = 35, hjust = 0, size = 11),
    axis.text.y       = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.box        = "vertical"
  )

ggsave(file.path(PLOT_DIR, "FigC_Evidence_Heatmap.pdf"),
       fig_C, width = 10, height = 8, dpi = 300)
ggsave(file.path(PLOT_DIR, "FigC_Evidence_Heatmap.png"),
       fig_C, width = 10, height = 8, dpi = 300)
message("Saved Figure C: Evidence integration heatmap")

# ── SUPPLEMENTARY TABLES ──────────────────────────────────────────────────────

## Supp Table: Full pMR results (all estimators)
supp_pMR <- pMR_all %>%
  left_join(select(pMR_ivw, protein, outcome, FDR, evidence_tier,
                   direction_consistent), by = c("protein","outcome")) %>%
  left_join(select(coloc, protein, PP.H4, coloc_call), by = "protein") %>%
  select(protein, outcome, method, n_snps, OR, CI_lo, CI_hi, pval,
         FDR, egger_intercept_p, cochran_Q_p,
         PP.H4, coloc_call, direction_consistent, evidence_tier) %>%
  arrange(protein, outcome, method)

write_csv(supp_pMR, file.path(OUT_DIR, "Supplementary_Table_S_pMR_FullResults.csv"))

## Supp Table: Coloc summary
supp_coloc <- coloc %>%
  left_join(fibroid_logfc, by = "protein") %>%
  select(protein, pQTL_id, region, n_snps,
         PP.H0, PP.H1, PP.H2, PP.H3, PP.H4,
         coloc_call, logfc_fib, direction)

write_csv(supp_coloc, file.path(OUT_DIR, "Supplementary_Table_S_Coloc_Summary.csv"))

message("\n=== All outputs saved ===")
message(sprintf("Figures : %s", PLOT_DIR))
message(sprintf("Tables  : %s", OUT_DIR))
message("\nFiles generated:")

message("  FigB_pMR_ForestPlot.pdf/png")
message("  FigC_Evidence_Heatmap.pdf/png")
message("  Supplementary_Table_S_pMR_FullResults.csv")
message("  Supplementary_Table_S_Coloc_Summary.csv")
