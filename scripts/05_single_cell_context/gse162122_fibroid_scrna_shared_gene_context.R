#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

root <- getwd()

data_dir <- file.path(root, "external_data", "GSE162122_fibroid_scrna")
out_dir <- file.path(root, "results", "07_single_cell_context")
asset_dir <- file.path(root, "results", "supplementary_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(asset_dir, showWarnings = FALSE, recursive = TRUE)

shared_path <- file.path(root, "data", "07_single_cell_context", "Supplementary_Table_S14c_Shared_Genes_MAGMA.tsv")
shared <- read.delim(shared_path, check.names = FALSE)
target_genes <- unique(shared$SYMBOL)

signatures <- list(
  Smooth_muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES", "MYLK", "TPM2"),
  Fibroblast_Stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "C7"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "CLDN5", "ESAM"),
  Pericyte = c("RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3"),
  Macrophage_Monocyte = c("LST1", "C1QA", "C1QB", "CD68", "CD14", "LYZ"),
  T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "KLRD1"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2")
)

find_matrix_sets <- function(base) {
  exact_mtx <- list.files(base, pattern = "^matrix[.]mtx[.]gz$", recursive = TRUE, full.names = TRUE)
  exact <- data.frame(
    Sample = vapply(dirname(exact_mtx), sample_label, character(1)),
    Matrix = exact_mtx,
    Features = file.path(dirname(exact_mtx), "features.tsv.gz"),
    stringsAsFactors = FALSE
  )

  top_mtx <- list.files(base, pattern = "_matrix[.]mtx[.]gz$", recursive = FALSE, full.names = TRUE)
  top_prefix <- sub("_matrix[.]mtx[.]gz$", "", top_mtx)
  top <- data.frame(
    Sample = basename(top_prefix),
    Matrix = top_mtx,
    Features = paste0(top_prefix, "_features.tsv.gz"),
    stringsAsFactors = FALSE
  )

  sets <- rbind(exact, top)
  sets <- sets[file.exists(sets$Matrix) & file.exists(sets$Features), ]
  sets[!duplicated(sets$Sample), ]
}

read_features <- function(f) {
  if (!file.exists(f)) stop("Missing features.tsv.gz: ", f)
  features <- read.delim(gzfile(f), header = FALSE, stringsAsFactors = FALSE)
  gene_symbols <- features[[2]]
  make.unique(gene_symbols)
}

sample_label <- function(path) {
  pieces <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1]]
  gsm <- pieces[grep("^GSM", pieces)]
  if (length(gsm) == 0) basename(path) else gsm[length(gsm)]
}

condition_label <- function(sample) {
  if (grepl("Fibroid", sample, ignore.case = TRUE)) "Fibroid" else "Myometrium"
}

score_celltypes <- function(log_norm, genes) {
  scores <- matrix(0, nrow = ncol(log_norm), ncol = length(signatures))
  colnames(scores) <- names(signatures)
  for (nm in names(signatures)) {
    idx <- which(genes %in% signatures[[nm]])
    if (length(idx) > 0) {
      scores[, nm] <- Matrix::colMeans(log_norm[idx, , drop = FALSE])
    }
  }
  best <- colnames(scores)[max.col(scores, ties.method = "first")]
  best_score <- apply(scores, 1, max)
  best[best_score <= 0] <- "Unassigned"
  best
}

summarise_sample <- function(row) {
  sample <- row$Sample
  condition <- condition_label(sample)
  message("Reading ", sample, " (", condition, ")")

  mtx <- Matrix::readMM(gzfile(row$Matrix))
  genes <- read_features(row$Features)
  rownames(mtx) <- genes

  totals <- Matrix::colSums(mtx)
  totals[totals == 0] <- 1
  log_norm <- log1p(t(t(mtx) / totals) * 10000)
  cell_type <- score_celltypes(log_norm, genes)

  present_targets <- intersect(target_genes, genes)
  rows <- list()
  i <- 1
  for (gene in target_genes) {
    if (!(gene %in% genes)) {
      rows[[i]] <- data.frame(
        Sample = sample, Condition = condition, Inferred_cell_type = "Not_in_matrix",
        Gene = gene, Cells = 0, Mean_log1pCP10K = NA_real_, Percent_expressing = NA_real_,
        stringsAsFactors = FALSE
      )
      i <- i + 1
      next
    }
    values <- as.numeric(log_norm[gene, ])
    raw <- as.numeric(mtx[gene, ])
    for (ct in sort(unique(cell_type))) {
      idx <- which(cell_type == ct)
      rows[[i]] <- data.frame(
        Sample = sample, Condition = condition, Inferred_cell_type = ct,
        Gene = gene, Cells = length(idx),
        Mean_log1pCP10K = mean(values[idx]),
        Percent_expressing = mean(raw[idx] > 0) * 100,
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
  counts <- as.data.frame(table(Inferred_cell_type = cell_type), stringsAsFactors = FALSE)
  counts$Sample <- sample
  counts$Condition <- condition
  counts <- counts[, c("Sample", "Condition", "Inferred_cell_type", "Freq")]
  colnames(counts)[4] <- "Cells"
  list(expr = do.call(rbind, rows), counts = counts)
}

matrix_sets <- find_matrix_sets(data_dir)
if (nrow(matrix_sets) == 0) stop("No matrix files found under ", data_dir)

results <- lapply(split(matrix_sets, seq_len(nrow(matrix_sets))), summarise_sample)
expr <- do.call(rbind, lapply(results, `[[`, "expr"))
counts <- do.call(rbind, lapply(results, `[[`, "counts"))

weighted <- aggregate(
  cbind(weighted_mean = expr$Mean_log1pCP10K * expr$Cells,
        weighted_pct = expr$Percent_expressing * expr$Cells,
        Cells = expr$Cells),
  by = list(Condition = expr$Condition, Inferred_cell_type = expr$Inferred_cell_type, Gene = expr$Gene),
  FUN = sum,
  na.rm = TRUE
)
weighted$Mean_log1pCP10K <- weighted$weighted_mean / weighted$Cells
weighted$Percent_expressing <- weighted$weighted_pct / weighted$Cells
weighted <- weighted[, c("Condition", "Inferred_cell_type", "Gene", "Cells", "Mean_log1pCP10K", "Percent_expressing")]
weighted <- merge(
  weighted,
  shared[, c("SYMBOL", "LOG10P_fib", "LOG10P_endo", "P_fib", "P_endo")],
  by.x = "Gene", by.y = "SYMBOL", all.x = TRUE
)
weighted <- weighted[order(weighted$Gene, weighted$Condition, -weighted$Mean_log1pCP10K), ]

write.csv(expr, file.path(out_dir, "GSE162122_Fibroid_scRNA_SharedGenes_BySample.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "GSE162122_Fibroid_scRNA_InferredCellType_Counts.csv"), row.names = FALSE)
write.csv(weighted, file.path(out_dir, "Table_GSE162122_Fibroid_scRNA_SharedGenes_CellContext.csv"), row.names = FALSE)
write.csv(weighted, file.path(asset_dir, "S14g_GSE162122_Fibroid_scRNA_SharedGenes_CellContext.csv"), row.names = FALSE)

message("Wrote GSE162122 fibroid scRNA shared-gene context tables")
