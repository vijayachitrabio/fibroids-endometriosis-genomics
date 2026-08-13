# ==============================================================================
# 02b_coloc_local.R
# Purpose : PARTIAL local version of Script 02.
#           The FIBROID GWAS side of coloc.abf is loaded from the local
#           FinnGen R9 file instead of the OpenGWAS API.
#           The pQTL side still requires one of:
#             (a) A valid OpenGWAS JWT (set OPENGWAS_JWT), OR
#             (b) Pre-downloaded regional pQTL .txt files in PQTL_DATA_DIR.
#
#           If OPENGWAS_JWT is set  → full coloc runs.
#           If not                  → pQTL regional extraction is skipped;
#                                     script produces a stub coloc_results
#                                     with PP.H4 = NA so downstream scripts
#                                     can still run (treated as "no coloc data").
#
# Outputs : coloc_results_summary.csv
#           coloc_results_all_posteriors.csv
#
# Author  : Vijayachitra Modhukur  |  April 2026
# Depends : coloc, ieugwasr, dplyr, readr, data.table
# ==============================================================================

library(coloc)
library(dplyr)
library(readr)
library(data.table)

# Optional: ieugwasr only loaded if token available
jwt_available <- nchar(Sys.getenv("OPENGWAS_JWT")) > 0
if (jwt_available) {
  library(ieugwasr)
  message("OPENGWAS_JWT found — pQTL extraction via API enabled.")
} else {
  message("OPENGWAS_JWT not set — running in fibroid-only mode.")
  message("PP.H4 will be NA unless pre-downloaded pQTL files are available.")
}

# ── PATHS ─────────────────────────────────────────────────────────────────────
BASE_DIR      <- Sys.getenv("FIBROID_BASE_DIR", unset = getwd())
OUT_DIR       <- Sys.getenv("FIBROID_OUT_DIR",  unset = file.path(BASE_DIR, "outputs_v2"))
GWAS_DIR      <- Sys.getenv("FINNGEN_GWAS_DIR",
                              unset = file.path(BASE_DIR, "gwas_scripts"))
PQTL_DATA_DIR <- Sys.getenv("PQTL_DATA_DIR",
                              unset = file.path(BASE_DIR, "pqtl_regional"))  # optional
dir.create(OUT_DIR, showWarnings = FALSE)

# ── FIBROID GWAS (local) ────────────────────────────────────────────────────
# Load the FinnGen R9 uterine myoma GWAS into memory once.
# ~20M rows; uses data.table for speed (~60 sec on most laptops).
FIBROID_FILE <- file.path(GWAS_DIR, "finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz")
message(sprintf("\nLoading local fibroid GWAS: %s", basename(FIBROID_FILE)))
message("(~20M rows, may take 30-90 sec) ...")

fibroid_gwas <- fread(FIBROID_FILE, sep = "\t", data.table = FALSE,
  col.names = c("chrom","pos","ref","alt","rsids","nearest_genes",
                "pval","mlogp","beta","sebeta","af_alt",
                "af_alt_cases","af_alt_controls"),
  skip = 1
)
fibroid_gwas$rsid_primary <- sub(",.*", "", fibroid_gwas$rsids)

message(sprintf("Loaded %d variants from fibroid GWAS.", nrow(fibroid_gwas)))

# ── pQTL IDs ──────────────────────────────────────────────────────────────────
pQTL_ids <- read_csv(file.path(OUT_DIR, "pQTL_GWAS_IDs.csv"),
                      show_col_types = FALSE) %>%
  filter(!is.na(opengwas_id))

instrument_summary <- read_csv(file.path(OUT_DIR, "pQTL_instrument_summary.csv"),
                                 show_col_types = FALSE)

# ── HELPER: extract fibroid GWAS region (from local file) ─────────────────────
pull_fibroid_region <- function(chrom, start_bp, end_bp) {
  subset(fibroid_gwas,
         chrom == as.integer(chrom) &
         pos   >= start_bp          &
         pos   <= end_bp)
}

# ── HELPER: extract pQTL region ───────────────────────────────────────────────
# Tries (1) pre-downloaded local file, then (2) OpenGWAS API if JWT available.
pull_pqtl_region <- function(gwas_id, chrom, start_bp, end_bp) {
  # Try local file first
  local_file <- file.path(PQTL_DATA_DIR,
                           sprintf("%s_chr%d_%d_%d.csv", gwas_id, chrom, start_bp, end_bp))
  if (file.exists(local_file)) {
    message("    [local] Reading pre-downloaded pQTL region.")
    return(read_csv(local_file, show_col_types = FALSE))
  }

  # Fall back to API
  if (!jwt_available) {
    message(sprintf("    [skip] No JWT and no local pQTL file for %s — coloc skipped.", gwas_id))
    return(NULL)
  }

  region_str <- sprintf("%d:%d-%d", chrom, start_bp, end_bp)
  tryCatch({
    res <- ieugwasr::associations(
      variants = region_str,
      id       = gwas_id,
      proxies  = 0
    )
    if (!is.null(res) && nrow(res) > 0) return(res)
    return(NULL)
  }, error = function(e) {
    message(sprintf("    [API] Region pull failed: %s", e$message))
    return(NULL)
  })
}

# ── HELPER: format coloc dataset ─────────────────────────────────────────────
format_coloc_dataset <- function(ss, beta_col, se_col, snp_col,
                                  type = "quant", N, sdY = 1) {
  ss <- ss[!is.na(ss[[beta_col]]) &
           !is.na(ss[[se_col]])   &
           ss[[se_col]] > 0       &
           !is.na(ss[[snp_col]]), ]
  ss <- ss[!duplicated(ss[[snp_col]]), ]
  if (nrow(ss) < 10) return(NULL)

  list(
    snp     = ss[[snp_col]],
    beta    = ss[[beta_col]],
    varbeta = ss[[se_col]]^2,
    type    = type,
    N       = N,
    sdY     = sdY
  )
}

# ── MAIN COLOC LOOP ───────────────────────────────────────────────────────────
N_fibroid <- 218000
N_decode  <- 35559
WINDOW_KB <- 500

coloc_results  <- list()
all_posteriors <- list()

for (i in seq_len(nrow(pQTL_ids))) {

  prot    <- pQTL_ids$protein[i]
  pqtl_id <- pQTL_ids$opengwas_id[i]

  lead_row <- instrument_summary %>% filter(protein == prot)
  if (nrow(lead_row) == 0) {
    message(sprintf("[%d/%d] No instrument info for %s — skipping.",
                    i, nrow(pQTL_ids), prot))
    next
  }

  chrom    <- as.integer(lead_row$lead_snp_chr[1])
  lead_pos <- as.integer(lead_row$lead_snp_pos[1])
  start_bp <- max(1, lead_pos - WINDOW_KB * 1000)
  end_bp   <- lead_pos + WINDOW_KB * 1000

  message(sprintf("\n[%d/%d] Coloc: %s | chr%d:%d-%d",
                  i, nrow(pQTL_ids), prot, chrom, start_bp, end_bp))

  # ── Fibroid GWAS region (local) ──
  ss_fibroid <- pull_fibroid_region(chrom, start_bp, end_bp)
  message(sprintf("  FinnGen fibroid region: %d variants", nrow(ss_fibroid)))

  if (nrow(ss_fibroid) < 50) {
    message(sprintf("  Too few fibroid variants (<50) in window. Skipping %s.", prot))
    next
  }

  # ── pQTL region (API or local) ──
  ss_pqtl <- pull_pqtl_region(pqtl_id, chrom, start_bp, end_bp)

  if (is.null(ss_pqtl) || nrow(ss_pqtl) < 50) {
    message(sprintf("  pQTL data unavailable for %s — storing NA coloc result.", prot))
    coloc_results[[prot]] <- data.frame(
      protein    = prot,
      pQTL_id    = pqtl_id,
      region     = sprintf("chr%d:%d-%d", chrom, start_bp, end_bp),
      n_snps     = NA_integer_,
      PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
      PP.H3 = NA_real_, PP.H4 = NA_real_,
      coloc_call = "Pending (pQTL data required)",
      stringsAsFactors = FALSE
    )
    next
  }

  # Harmonise on shared rsids
  fib_rsid  <- if ("rsid_primary" %in% names(ss_fibroid)) ss_fibroid$rsid_primary
               else sub(",.*", "", ss_fibroid$rsids)
  pqtl_rsid <- if ("rsid" %in% names(ss_pqtl)) ss_pqtl$rsid
               else ss_pqtl$SNP

  common <- intersect(fib_rsid, pqtl_rsid)
  message(sprintf("  Shared SNPs: %d", length(common)))
  if (length(common) < 50) {
    message("  Too few shared SNPs (<50). Skipping.")
    next
  }

  ss_fib_f <- ss_fibroid[fib_rsid %in% common, ]
  ss_fib_f$SNP <- sub(",.*", "", ss_fib_f$rsids)

  ss_pqtl_f <- ss_pqtl[pqtl_rsid %in% common, ]
  if ("rsid" %in% names(ss_pqtl_f))  ss_pqtl_f$SNP <- ss_pqtl_f$rsid
  if ("beta" %in% names(ss_pqtl_f))  ss_pqtl_f$BETA <- ss_pqtl_f$beta
  if ("se"   %in% names(ss_pqtl_f))  ss_pqtl_f$SE   <- ss_pqtl_f$se

  # Format for coloc
  d1 <- format_coloc_dataset(ss_pqtl_f, "BETA", "SE", "SNP",
                               type = "quant", N = N_decode, sdY = 1)
  d2 <- format_coloc_dataset(ss_fib_f,  "beta", "sebeta", "SNP",
                               type = "cc",    N = N_fibroid)

  if (is.null(d1) || is.null(d2)) {
    message("  Could not format coloc datasets. Skipping.")
    next
  }

  # Run coloc.abf
  tryCatch({
    result <- coloc.abf(dataset1 = d1, dataset2 = d2,
                         p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

    summary_row <- data.frame(
      protein    = prot,
      pQTL_id    = pqtl_id,
      region     = sprintf("chr%d:%d-%d", chrom, start_bp, end_bp),
      n_snps     = result$summary["nsnps"],
      PP.H0      = round(result$summary["PP.H0.abf"], 4),
      PP.H1      = round(result$summary["PP.H1.abf"], 4),
      PP.H2      = round(result$summary["PP.H2.abf"], 4),
      PP.H3      = round(result$summary["PP.H3.abf"], 4),
      PP.H4      = round(result$summary["PP.H4.abf"], 4),
      coloc_call = case_when(
        result$summary["PP.H4.abf"] >= 0.80 ~ "Strong (H4≥0.80)",
        result$summary["PP.H4.abf"] >= 0.50 ~ "Moderate (H4 0.50-0.79)",
        result$summary["PP.H4.abf"] >= 0.20 ~ "Weak (H4 0.20-0.49)",
        TRUE                                 ~ "No evidence (H4<0.20)"
      ),
      stringsAsFactors = FALSE
    )
    coloc_results[[prot]] <- summary_row

    posterior_df <- result$results %>%
      mutate(protein = prot,
             region  = sprintf("chr%d:%d-%d", chrom, start_bp, end_bp))
    all_posteriors[[prot]] <- posterior_df

    message(sprintf("  PP.H4 = %.3f (%s)",
                    result$summary["PP.H4.abf"], summary_row$coloc_call))

  }, error = function(e) {
    message(sprintf("  coloc ERROR for %s: %s", prot, e$message))
  })
}

# ── SAVE ──────────────────────────────────────────────────────────────────────
coloc_summary <- bind_rows(coloc_results) %>% arrange(desc(PP.H4))

write_csv(coloc_summary, file.path(OUT_DIR, "coloc_results_summary.csv"))

if (length(all_posteriors) > 0) {
  write_csv(bind_rows(all_posteriors),
            file.path(OUT_DIR, "coloc_results_all_posteriors.csv"))
}

message("\n=== Colocalization Summary ===")
print(select(coloc_summary, protein, PP.H4, coloc_call))
message(sprintf("\nProteins with PP.H4 data: %d",  sum(!is.na(coloc_summary$PP.H4))))
message(sprintf("Proteins with PP.H4 ≥ 0.5 : %d", sum(coloc_summary$PP.H4 >= 0.5, na.rm=TRUE)))
message("\nColoc analysis complete. Output: coloc_results_summary.csv")
message()
message("If PP.H4 = NA for all proteins: OpenGWAS JWT is needed for pQTL data.")
message("Renew token at: https://api.opengwas.io/account")
message("Then set:  echo 'OPENGWAS_JWT=eyJ...' >> ~/.Renviron  and restart R.")
