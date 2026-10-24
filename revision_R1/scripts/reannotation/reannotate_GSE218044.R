# Cell-level reference-based reannotation of GSE218044 (scEndoExplorer adenomyosis dataset)
# and re-derivation of S-LDSC marker genes from confidently labelled cells.
#
# 1. Per-cell SingleR (HPCA and Blueprint, label.main), same reference choice and
#    label harmonisation as scEndoExplorer's cluster-level annotation.
# 2. Per-cell confidence: SingleR pruned label (outlier on delta) and delta.next margin;
#    concordance of the per-cell label with the original cluster-level label.
# 3. Markers with scEndoExplorer's one_vs_rest_wilcox (precalculate.R):
#      (a) all cells, original labels  -> must reproduce cell_type_marker_genes.csv
#      (b) confident cells only (per-cell label == cluster label, not pruned)
#    then the same selection as run_ldsc_seg.sh: log2FC > 0.25, padj < 0.05,
#    top 500 per cell type by padj then log2FC.
suppressPackageStartupMessages({
  library(SingleR); library(celldex); library(Matrix); library(BiocParallel)
  library(data.table); library(dplyr)
})
SC <- "/Users/vijayachitramodhukur/Desktop/ECLAI/scEndoExplorer_v1"
OUT <- getwd()
t0 <- Sys.time(); tick <- function(...) cat(format(Sys.time() - t0), "|", ..., "\n")

x <- readRDS(file.path(SC, "data/adenomyosis_endo.rds"))
counts <- x$expr
meta <- as.data.table(x$metadata)[match(colnames(counts), cell_id)]
stopifnot(all(meta$cell_id == colnames(counts)))
tick("cells", ncol(counts), "genes", nrow(counts))

# label map and Wilcoxon function, taken verbatim from scEndoExplorer
src <- readLines(file.path(SC, "data_prep/annotate_adenomyosis.R"))
i0 <- grep("^UTERUS_LABEL_MAP <- c\\(", src); i1 <- grep("^harmonise <- function", src)
eval(parse(text = src[i0:(i1 - 1)]))
pc <- readLines(file.path(SC, "data_prep/precalculate.R"))
j0 <- grep("^one_vs_rest_wilcox <- function", pc); j1 <- grep("^# ---- 3. Pre-compute", pc)
eval(parse(text = pc[j0:(j1 - 1)]))
harm <- function(l) { o <- unname(UTERUS_LABEL_MAP[l]); o[is.na(o)] <- "Other"; o }

# ---------------------------------------------------------------- 1. per-cell SingleR
lib <- Matrix::colSums(counts)
lognorm <- counts %*% Diagonal(x = 1e4 / lib)
lognorm@x <- log1p(lognorm@x)
dimnames(lognorm) <- dimnames(counts)
hpca <- HumanPrimaryCellAtlasData(); bp <- BlueprintEncodeData()
bpp <- MulticoreParam(4)
tick("SingleR HPCA ...")
rh <- SingleR(lognorm, hpca, hpca$label.main, BPPARAM = bpp)
tick("SingleR Blueprint ...")
rb <- SingleR(lognorm, bp, bp$label.main, BPPARAM = bpp)
sh <- apply(rh$scores, 1, max); sb <- apply(rb$scores, 1, max)
use_b <- sb >= sh
cell <- data.table(
  cell_id = meta$cell_id, donor = meta$patient_id, cluster = meta$cluster, cluster_label = meta$cell_type,
  ref = ifelse(use_b, "Blueprint", "HPCA"),
  raw_label = ifelse(use_b, rb$labels, rh$labels),
  pruned = ifelse(use_b, is.na(rb$pruned.labels), is.na(rh$pruned.labels)),
  delta_next = ifelse(use_b, rb$delta.next, rh$delta.next),
  max_score = pmax(sh, sb))
cell[, cell_label := harm(raw_label)]
cell[, status := fifelse(pruned, "low confidence (pruned)",
                  fifelse(cell_label == cluster_label, "concordant", "discordant"))]
fwrite(cell, file.path(OUT, "GSE218044_cell_level_labels.csv.gz"))
tick("labels written")

summ_ct <- cell[, .(cells = .N, concordant = mean(status == "concordant"),
                    discordant = mean(status == "discordant"), pruned = mean(pruned),
                    median_delta = median(delta_next)), by = cluster_label][order(cluster_label)]
summ_donor <- cell[, .(cells = .N, concordant = mean(status == "concordant"),
                       discordant = mean(status == "discordant"), pruned = mean(pruned)), by = .(donor, cluster_label)]
conf_mat <- dcast(cell[, .N, by = .(cluster_label, cell_label)], cluster_label ~ cell_label, value.var = "N", fill = 0)
fwrite(summ_ct, file.path(OUT, "reannotation_summary_by_celltype.csv"))
fwrite(summ_donor, file.path(OUT, "reannotation_summary_by_donor.csv"))
fwrite(conf_mat, file.path(OUT, "reannotation_confusion_cluster_vs_cell.csv"))
print(summ_ct)

# ---------------------------------------------------------------- 2-3. markers
select_top <- function(mk) {
  mk <- mk[!is.na(mk$gene) & !is.na(mk$avg_log2FC) & !is.na(mk$p_val_adj), ]
  mk <- mk[mk$avg_log2FC > 0.25 & mk$p_val_adj < 0.05, ]
  as.data.frame(mk %>% group_by(cell_type) %>% arrange(p_val_adj, desc(avg_log2FC)) %>%
                  slice_head(n = 500) %>% ungroup())[, c("gene", "cell_type", "avg_log2FC", "p_val_adj")]
}
markers_for <- function(labels, cells) {
  ids <- colnames(counts)
  cts <- sort(unique(na.omit(labels)))
  do.call(rbind, lapply(cts, function(ct) {
    fg <- intersect(cells[labels == ct], ids); bg <- setdiff(intersect(cells, ids), fg)
    df <- one_vs_rest_wilcox(counts, fg, bg)
    if (nrow(df)) { df$cluster <- ct; df$cell_type <- ct }
    tick("  markers", ct, nrow(df)); df
  }))
}
orig <- fread(file.path(SC, "data/cell_type_marker_genes.csv"))
tick("markers (a): all cells, original labels")
ma <- select_top(markers_for(cell$cluster_label, cell$cell_id))
fwrite(ma, file.path(OUT, "markers_all_cells_reproduction.csv"))
rep_ok <- nrow(ma) == nrow(orig) && all(ma$gene == orig$gene) && all(ma$cell_type == orig$cell_type)
tick("reproduction of cell_type_marker_genes.csv exact:", rep_ok,
     "| per-type overlap:", paste(sapply(unique(orig$cell_type), function(ct) length(intersect(ma$gene[ma$cell_type == ct], orig$gene[orig$cell_type == ct]))), collapse = ","))
conf <- cell[status == "concordant"]
tick("markers (b): confident cells only,", nrow(conf), "cells")
mb <- select_top(markers_for(conf$cluster_label, conf$cell_id))
fwrite(mb, file.path(OUT, "markers_confident_cells.csv"))
ov <- rbindlist(lapply(sort(unique(orig$cell_type)), function(ct) {
  a <- orig$gene[orig$cell_type == ct]; b <- mb$gene[mb$cell_type == ct]
  data.table(cell_type = ct, n_original = length(a), n_confident = length(b), shared = length(intersect(a, b)),
             jaccard = round(length(intersect(a, b)) / length(union(a, b)), 3))
}))
fwrite(ov, file.path(OUT, "marker_overlap_original_vs_confident.csv")); print(ov)
shared10 <- c("SYNE1","GREB1","WNT4","ESR1","EFR3B","VIP","LDLRAD2","DNM3","PLEKHM3","ARL14EP")
hit <- function(m) sapply(shared10, function(g) paste(m$cell_type[m$gene == g], collapse = ";"))
fwrite(data.table(gene = shared10, original = hit(orig), confident_cells = hit(mb)), file.path(OUT, "shared_gene_marker_status.csv"))
tick("done")
