library(ggplot2)
library(dplyr)
library(stringr)

# Get repo root correctly when sourced or Rscript
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  BASE_DIR <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) getwd())
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    BASE_DIR <- dirname(normalizePath(sub("--file=", "", file_arg[1])))
  } else {
    BASE_DIR <- getwd()
  }
}
repo_root <- normalizePath(file.path(BASE_DIR, "../.."), mustWork = FALSE)

# Load data
res_file <- file.path(repo_root, "results/PheWAS_MR_All_Results.csv")
df <- read.csv(res_file, stringsAsFactors = FALSE)

# Filter for FDR significant and exclude cancers
cancers <- c("endometrial_cancer", "breast_cancer", "ovarian_cancer", "cervical_cancer", "leiomyosarcoma")
plot_data <- df %>%
  filter(Conclusion == "FDR-significant") %>%
  filter(!Outcome %in% cancers)

# Format names
format_name <- function(x) {
  x <- str_replace_all(x, "_", " ")
  x <- tools::toTitleCase(x)
  # Fix specific acronyms
  x <- str_replace(x, "Pcos", "PCOS")
  return(x)
}

plot_data <- plot_data %>%
  mutate(
    Outcome_clean = format_name(Outcome),
    Exposure_clean = ifelse(Exposure == "Uterine_Fibroids", "Uterine fibroids", "Endometriosis"),
    Label = paste(Exposure_clean, "->", Outcome_clean),
    # Parse Confidence Interval bounds
    CI_lower = as.numeric(sapply(strsplit(IVW_CI, "–"), `[`, 1)),
    CI_upper = as.numeric(sapply(strsplit(IVW_CI, "–"), `[`, 2))
  )

# Sort by Exposure then OR
plot_data <- plot_data %>%
  arrange(Exposure_clean, desc(IVW_OR)) %>%
  mutate(Label = factor(Label, levels = rev(Label)))

# Add N_instruments to legend labels
exp_n <- plot_data %>%
  group_by(Exposure_clean) %>%
  summarize(N = first(N_instruments))

legend_labels <- c(
  "Endometriosis" = paste0("Endometriosis instruments (n=", exp_n$N[exp_n$Exposure_clean == "Endometriosis"], ")"),
  "Uterine fibroids" = paste0("Uterine fibroids instruments (n=", exp_n$N[exp_n$Exposure_clean == "Uterine fibroids"], ")")
)

# Colors and shapes matching image
colors <- c("Uterine fibroids" = "#00897B", "Endometriosis" = "#E65100")
shapes <- c("Uterine fibroids" = 15, "Endometriosis" = 16)

# Create Plot
p <- ggplot(plot_data, aes(x = IVW_OR, y = Label, color = Exposure_clean, shape = Exposure_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.15, size = 1) +
  geom_point(size = 4.5) +
  geom_text(aes(label = sprintf("OR=%.3f (%.3f–%.3f); FDR<%.3f", 
                                IVW_OR, CI_lower, CI_upper, 
                                ifelse(IVW_FDR_P < 0.001, 0.001, IVW_FDR_P)),
                x = CI_upper + 0.03), 
            hjust = 0, size = 3, color = "gray30", show.legend = FALSE) +
  scale_color_manual(values = colors, labels = legend_labels) +
  scale_shape_manual(values = shapes, labels = legend_labels) +
  scale_x_continuous(limits = c(0.7, max(plot_data$CI_upper) + 0.5)) +
  labs(
    x = "Odds ratio (OR) for exposure -> outcome, with 95% confidence interval",
    y = NULL,
    color = "MR exposure",
    shape = "MR exposure"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold")
  )

# Save
out_dir <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "Figure_PheWAS_MR_Summary_NoCancers.png"), p, width = 11, height = 3.5, dpi = 600, bg="white")
ggsave(file.path(out_dir, "Figure_PheWAS_MR_Summary_NoCancers.pdf"), p, width = 11, height = 3.5, dpi = 600, bg="white")

cat("Saved high-res figures to", out_dir, "\n")
