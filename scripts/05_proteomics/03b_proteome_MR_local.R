# ==============================================================================
# 03b_proteome_MR_local.R
# Purpose : NO-API version of Script 03.
#           Identical analysis (proteome-wide MR, 16 proteins × 7 outcomes),
#           but outcome data is extracted from LOCAL FinnGen R9 summary stats
#           files instead of the OpenGWAS API.
#
#           API dependency removed:
#             extract_outcome_data() → replaced by extract_outcome_local()
#             which reads from pre-downloaded FinnGen .gz files.
#
#           Coeliac disease (originally ebi-a-GCST005523, UKBB) is dropped
#           from this local run because no FinnGen equivalent is available.
#           Replace with finn-b-K11_COELIAC when/if token is restored.
#
#           Key inputs:
#             pQTL_instruments_all.csv  — from Script 01 or 01b
#             pQTL_GWAS_IDs.csv         — protein → OpenGWAS ID mapping
#             coloc_results_summary.csv — from Script 02 (or 02b)
#             finngen_R9_*.gz           — local FinnGen R9 summary stats
#
# Outputs : pMR_all_results.csv
#           pMR_ivw_results.csv
#           pMR_FDR_significant.csv
#           pMR_evidence_matrix.csv
#
# Author  : Vijayachitra Modhukur  |  April 2026
# Depends : TwoSampleMR, dplyr, readr, tidyr, data.table
# ==============================================================================

library(TwoSampleMR)
library(dplyr)
library(readr)
library(tidyr)
library(data.table)

# ── PATHS ─────────────────────────────────────────────────────────────────────
BASE_DIR <- Sys.getenv("FIBROID_BASE_DIR", unset = getwd())
OUT_DIR  <- Sys.getenv("FIBROID_OUT_DIR",  unset = file.path(BASE_DIR, "outputs_v2"))
GWAS_DIR <- Sys.getenv("FINNGEN_GWAS_DIR",
                        unset = file.path(BASE_DIR, "gwas_scripts"))

message("=== Script 03b: Proteome-MR (local FinnGen mode) ===")
message(sprintf("OUT_DIR : %s", OUT_DIR))
message(sprintf("GWAS_DIR: %s", GWAS_DIR))

# ── OUTCOMES ──────────────────────────────────────────────────────────────────
# 6/7 outcomes available locally (coeliac dropped — no FinnGen R9 equivalent)
outcomes <- data.frame(
  outcome_label  = c("PCOS",
                     "Ovarian cyst",
                     "Ovarian cancer",
                     "Hypertension",
                     "Depression",
                     "Breast cancer"),
  finngen_file   = c("finngen_R9_E4_PCOS.gz",
                     "finngen_R9_N14_OVARYCYST.gz",
                     "finngen_R9_C3_OVARY_EXALLC.gz",
                     "finngen_R9_I9_HYPTENS.gz",
                     "finngen_R9_F5_DEPRESSIO.gz",
                     "finngen_R9_C3_BREAST_EXALLC.gz"),
  outcome_N      = c(218882, 218882, 218882, 218882, 218882, 218882),
  stringsAsFactors = FALSE
)
message(sprintf("Running with %d outcomes (coeliac excluded — no local file)",
                nrow(outcomes)))

# ── INSTRUMENTS ───────────────────────────────────────────────────────────────
instruments_all <- read_csv(file.path(OUT_DIR, "pQTL_instruments_all.csv"),
                             show_col_types = FALSE)
pQTL_ids        <- read_csv(file.path(OUT_DIR, "pQTL_GWAS_IDs.csv"),
                             show_col_types = FALSE) %>%
  filter(!is.na(opengwas_id) | !is.na(protein))

# Colocalization results (may be absent if Script 02/02b not yet run)
coloc_path <- file.path(OUT_DIR, "coloc_results_summary.csv")
if (file.exists(coloc_path)) {
  coloc_summary <- read_csv(coloc_path, show_col_types = FALSE)
  message("Coloc results loaded.")
} else {
  warning("coloc_results_summary.csv not found — PP.H4 will be NA. ",
          "Run Script 02 or 02b first for full evidence tiering.")
  coloc_summary <- data.frame(protein = character(0),
                               PP.H4 = numeric(0),
                               coloc_call = character(0))
}

# Olink logFC direction
fibroid_logfc <- data.frame(
  protein   = c("CHRDL2","TNFRSF11B","EFEMP1","FLT3LG","EDA2R",
                "CDH3","COL9A1","FGF23","TNFRSF17","PTK7",
                "CD55","TNFSF11","LEFTY2","LPL","CHL1","TFPI"),
  logfc_fib = c(0.148,-0.059,-0.052,-0.053,-0.054,
                0.062,-0.074, 0.118,-0.058, 0.055,
               -0.033, 0.082, 0.078,-0.057,-0.036,-0.037),
  stringsAsFactors = FALSE
)

# ── LOCAL OUTCOME EXTRACTION ───────────────────────────────────────────────────
# Reads a FinnGen R9 .gz file and extracts rows matching the requested SNP IDs.
# FinnGen column format:
#   #chrom  pos  ref  alt  rsids  nearest_genes  pval  mlogp  beta  sebeta  af_alt  af_alt_cases  af_alt_controls
#
# Strategy: read the full file with data.table fread (fast), then match by:
#   1) exact rsID token in FinnGen's comma-separated rsids field
#   2) chromosome:position fallback, for lifted coordinate matches
# For files with ~20M rows this takes ~30-60 seconds per file.
# Cached per outcome across all protein runs within the session.

.outcome_cache <- new.env(parent = emptyenv())
.extraction_qc <- list()

regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

first_rsid_token <- function(rsids) {
  sub(",.*", "", rsids)
}

complement_allele <- function(x) {
  recode(x, A = "T", T = "A", C = "G", G = "C", .default = x)
}

extract_outcome_local <- function(instruments_df, outcome_label, finngen_file, outcome_N) {

  cache_key <- outcome_label

  # Load and cache full outcome file (once per outcome per session)
  if (!exists(cache_key, envir = .outcome_cache)) {
    fpath <- file.path(GWAS_DIR, finngen_file)
    if (!file.exists(fpath)) {
      warning(sprintf("FinnGen file not found: %s", fpath))
      return(NULL)
    }
    message(sprintf("  [Cache] Loading %s ...", finngen_file))
    dt <- fread(fpath, sep = "\t", data.table = FALSE,
                col.names = c("chrom","pos","ref","alt","rsids",
                               "nearest_genes","pval","mlogp","beta",
                               "sebeta","af_alt","af_alt_cases","af_alt_controls"),
                skip = 1)
    # Keep all FinnGen rsID aliases, but record a primary token for reporting.
    dt$rsid_primary <- first_rsid_token(dt$rsids)
    dt$coord_key <- paste(dt$chrom, dt$pos, sep = ":")
    assign(cache_key, dt, envir = .outcome_cache)
    message(sprintf("  [Cache] %s loaded (%d rows)", outcome_label, nrow(dt)))
  }

  dt <- get(cache_key, envir = .outcome_cache)

  snps <- instruments_df$SNP
  coord_keys <- paste(instruments_df$chr.exposure, instruments_df$pos.exposure, sep = ":")

  rsid_pattern <- paste0("(^|,)(", paste(regex_escape(snps), collapse = "|"), ")(,|$)")
  rsid_hit <- grepl(rsid_pattern, dt$rsids)
  coord_hit <- dt$coord_key %in% coord_keys

  candidate_rows <- dt[rsid_hit | coord_hit, ]

  if (nrow(candidate_rows) == 0) {
    message(sprintf("    No matching SNPs found in %s", outcome_label))
    .extraction_qc[[paste(outcome_label, paste(snps, collapse = "_"), sep = "::")]] <<-
      data.frame(
        outcome = outcome_label,
        protein = unique(instruments_df$protein),
        SNP = snps,
        requested_chr = instruments_df$chr.exposure,
        requested_pos = instruments_df$pos.exposure,
        match_status = "not_found",
        match_method = NA_character_,
        finngen_rsid = NA_character_,
        finngen_chr = NA_integer_,
        finngen_pos = NA_integer_,
        finngen_ref = NA_character_,
        finngen_alt = NA_character_,
        allele_match_class = NA_character_,
        stringsAsFactors = FALSE
      )
    return(NULL)
  }

  matched_list <- lapply(seq_len(nrow(instruments_df)), function(i) {
    inst <- instruments_df[i, ]
    by_rsid <- grepl(paste0("(^|,)", regex_escape(inst$SNP), "(,|$)"),
                     candidate_rows$rsids)
    by_coord <- candidate_rows$coord_key ==
      paste(inst$chr.exposure, inst$pos.exposure, sep = ":")

    rows <- candidate_rows[by_rsid | by_coord, ]
    if (nrow(rows) == 0) return(NULL)

    rows$SNP <- inst$SNP
    rows$requested_chr <- inst$chr.exposure
    rows$requested_pos <- inst$pos.exposure
    rows$match_method <- ifelse(by_rsid[by_rsid | by_coord], "rsid", "coordinate")
    rows$match_method[by_rsid[by_rsid | by_coord] & by_coord[by_rsid | by_coord]] <-
      "rsid+coordinate"
    rows
  })

  out_rows <- bind_rows(matched_list) %>%
    group_by(SNP) %>%
    arrange(match(match_method, c("rsid+coordinate", "rsid", "coordinate"))) %>%
    slice(1) %>%
    ungroup()

  allele_qc <- instruments_df %>%
    select(protein, SNP, chr.exposure, pos.exposure,
           effect_allele.exposure, other_allele.exposure) %>%
    left_join(
      out_rows %>%
        transmute(SNP,
                  match_method,
                  finngen_rsid = rsids,
                  finngen_chr = chrom,
                  finngen_pos = pos,
                  finngen_ref = ref,
                  finngen_alt = alt),
      by = "SNP"
    ) %>%
    mutate(
      match_status = if_else(is.na(finngen_rsid), "not_found", "found"),
      effect_allele.exposure_comp = complement_allele(effect_allele.exposure),
      other_allele.exposure_comp = complement_allele(other_allele.exposure),
      allele_match_class = case_when(
        is.na(finngen_rsid) ~ NA_character_,
        effect_allele.exposure == finngen_alt &
          other_allele.exposure == finngen_ref ~ "aligned",
        effect_allele.exposure == finngen_ref &
          other_allele.exposure == finngen_alt ~ "reversed",
        effect_allele.exposure_comp == finngen_alt &
          other_allele.exposure_comp == finngen_ref ~ "strand_aligned",
        effect_allele.exposure_comp == finngen_ref &
          other_allele.exposure_comp == finngen_alt ~ "strand_reversed",
        effect_allele.exposure %in% c("A", "T") &
          other_allele.exposure %in% c("A", "T") &
          finngen_alt %in% c("A", "T") &
          finngen_ref %in% c("A", "T") ~ "palindromic_or_ambiguous",
        effect_allele.exposure %in% c("C", "G") &
          other_allele.exposure %in% c("C", "G") &
          finngen_alt %in% c("C", "G") &
          finngen_ref %in% c("C", "G") ~ "palindromic_or_ambiguous",
        TRUE ~ "incompatible"
      ),
      outcome = outcome_label,
      requested_chr = chr.exposure,
      requested_pos = pos.exposure
    ) %>%
    select(outcome, protein, SNP, requested_chr, requested_pos,
           effect_allele.exposure, other_allele.exposure,
           match_status, match_method, finngen_rsid, finngen_chr, finngen_pos,
           finngen_ref, finngen_alt, allele_match_class)

  .extraction_qc[[paste(outcome_label, paste(snps, collapse = "_"), sep = "::")]] <<-
    allele_qc

  # Format to TwoSampleMR outcome format
  out_dat <- data.frame(
    SNP                   = out_rows$SNP,
    beta.outcome          = out_rows$beta,
    se.outcome            = out_rows$sebeta,
    effect_allele.outcome = out_rows$alt,    # FinnGen: alt = effect allele
    other_allele.outcome  = out_rows$ref,
    eaf.outcome           = out_rows$af_alt,
    pval.outcome          = out_rows$pval,
    samplesize.outcome    = outcome_N,
    outcome               = outcome_label,
    id.outcome            = outcome_label,
    mr_keep.outcome       = TRUE,
    stringsAsFactors      = FALSE
  )
  return(out_dat)
}

# ── MR FUNCTION ───────────────────────────────────────────────────────────────
run_proteome_mr_local <- function(protein_name, outcome_label,
                                  finngen_file, outcome_N,
                                  instruments_df) {

  # Filter instruments for this protein
  exp_dat <- instruments_df %>%
    filter(protein == protein_name)

  if (nrow(exp_dat) == 0) {
    message(sprintf("    No instruments for %s", protein_name))
    return(NULL)
  }

  # Convert to TwoSampleMR exposure format
  exp_dat <- exp_dat %>%
    mutate(
      exposure             = protein_name,
      id.exposure          = paste0("prot_", protein_name),
      units.exposure       = "SD",
      phenotype.col        = "exposure",
      mr_keep.exposure     = TRUE,
      samplesize.exposure  = 35559
    )

  # Extract outcome data from local file
  out_dat <- extract_outcome_local(exp_dat, outcome_label, finngen_file, outcome_N)
  if (is.null(out_dat) || nrow(out_dat) == 0) return(NULL)

  # Harmonise
  tryCatch({
    dat <- harmonise_data(exp_dat, out_dat, action = 2)
  }, error = function(e) {
    message(sprintf("    Harmonise error for %s→%s: %s", protein_name, outcome_label, e$message))
    return(NULL)
  })

  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat <- dat %>% filter(mr_keep)
  if (nrow(dat) == 0) return(NULL)

  n_snps <- nrow(dat)

  # Choose estimators based on instrument count
  methods <- if (n_snps >= 3) {
    c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression")
  } else if (n_snps == 2) {
    c("mr_ivw_mre", "mr_weighted_median")
  } else {
    c("mr_wald_ratio")   # single instrument → Wald ratio
  }

  mr_res <- tryCatch(
    mr(dat, method_list = methods),
    error = function(e) { message("    MR error: ", e$message); NULL }
  )
  if (is.null(mr_res) || nrow(mr_res) == 0) return(NULL)

  # Pleiotropy / heterogeneity (only meaningful with ≥3 instruments)
  egger_int_p <- NA_real_
  cochran_q_p <- NA_real_
  if (n_snps >= 3) {
    pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    het   <- tryCatch(mr_heterogeneity(dat),   error = function(e) NULL)
    if (!is.null(pleio) && nrow(pleio) > 0) egger_int_p <- pleio$pval[1]
    if (!is.null(het)   && nrow(het) > 0) {
      q_row <- het[grepl("^Inverse variance weighted", het$method), ]
      if (nrow(q_row) > 0) cochran_q_p <- q_row$Q_pval[1]
    }
  }

  # Format row per method
  fmt_row <- function(method_str, label) {
    if (method_str == "Inverse variance weighted") {
      row <- mr_res %>% filter(grepl("^Inverse variance weighted", method))
    } else {
      row <- mr_res %>% filter(method == method_str)
    }
    if (nrow(row) == 0) return(data.frame())
    data.frame(
      protein             = protein_name,
      outcome             = outcome_label,
      method              = label,
      n_snps              = n_snps,
      OR                  = round(exp(row$b),   4),
      CI_lo               = round(exp(row$b - 1.96 * row$se), 4),
      CI_hi               = round(exp(row$b + 1.96 * row$se), 4),
      beta                = round(row$b,         5),
      se                  = round(row$se,        5),
      pval                = signif(row$pval,     4),
      egger_intercept_p   = round(egger_int_p,  4),
      cochran_Q_p         = round(cochran_q_p,  4),
      stringsAsFactors    = FALSE
    )
  }

  bind_rows(
    fmt_row("Inverse variance weighted",  "IVW"),
    fmt_row("Wald ratio",                 "IVW"),   # single-SNP
    fmt_row("Weighted median",            "WM"),
    fmt_row("MR Egger",                  "Egger")
  )
}

# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
proteins_to_run <- unique(instruments_all$protein)
message(sprintf("\nRunning proteome-MR: %d proteins × %d outcomes = %d tests",
                length(proteins_to_run), nrow(outcomes),
                length(proteins_to_run) * nrow(outcomes)))

results_list <- list()

for (prot in proteins_to_run) {
  message(sprintf("\n--- Protein: %s ---", prot))

  for (j in seq_len(nrow(outcomes))) {
    out_label <- outcomes$outcome_label[j]
    out_file  <- outcomes$finngen_file[j]
    out_N     <- outcomes$outcome_N[j]
    message(sprintf("  → %s", out_label))

    res <- run_proteome_mr_local(
      protein_name  = prot,
      outcome_label = out_label,
      finngen_file  = out_file,
      outcome_N     = out_N,
      instruments_df = instruments_all
    )
    if (!is.null(res) && nrow(res) > 0) {
      results_list[[paste0(prot, "_", out_label)]] <- res
    }
  }
}

# ── COMPILE AND ANNOTATE ──────────────────────────────────────────────────────
pMR_all <- bind_rows(results_list)

if (nrow(pMR_all) == 0) {
  stop("No MR results produced. Check instruments file and FinnGen GWAS paths.")
}

# FDR on IVW results, within each outcome
pMR_ivw <- pMR_all %>%
  filter(method == "IVW") %>%
  group_by(outcome) %>%
  mutate(FDR = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  left_join(fibroid_logfc, by = "protein") %>%
  left_join(
    if (nrow(coloc_summary) > 0) select(coloc_summary, protein, PP.H4, coloc_call)
    else data.frame(protein=character(), PP.H4=numeric(), coloc_call=character()),
    by = "protein"
  ) %>%
  mutate(
    direction_consistent = case_when(
      logfc_fib > 0 & beta > 0 ~ "Consistent (both positive)",
      logfc_fib < 0 & beta < 0 ~ "Consistent (both negative)",
      TRUE                      ~ "Discordant"
    ),
    evidence_tier = case_when(
      FDR < 0.05 & !is.na(PP.H4) & PP.H4 >= 0.80 & direction_consistent != "Discordant" ~
        "Tier 1: FDR-sig + Coloc + Consistent",
      FDR < 0.05 & !is.na(PP.H4) & PP.H4 >= 0.50 & direction_consistent != "Discordant" ~
        "Tier 2: FDR-sig + Mod.Coloc + Consistent",
      FDR < 0.05 & direction_consistent != "Discordant" ~
        "Tier 3: FDR-sig + Consistent (no/pending coloc)",
      FDR < 0.05 ~
        "Tier 4: FDR-sig only (direction discordant or no coloc)",
      TRUE ~ "Non-significant"
    )
  ) %>%
  arrange(FDR, pval)

pMR_sig <- pMR_ivw %>% filter(FDR < 0.05)
message(sprintf("\nFDR-significant proteome-MR hits: %d", nrow(pMR_sig)))
if (nrow(pMR_sig) > 0) {
  print(select(pMR_sig, protein, outcome, OR, CI_lo, CI_hi, pval, FDR, evidence_tier))
}

# Evidence matrix (wide format for heatmap)
evidence_matrix <- pMR_ivw %>%
  mutate(OR_display = ifelse(FDR < 0.05, OR, NA)) %>%
  select(protein, outcome, OR_display, FDR, PP.H4, evidence_tier) %>%
  pivot_wider(names_from = outcome,
              values_from = c(OR_display, FDR, PP.H4, evidence_tier))

# ── SAVE ──────────────────────────────────────────────────────────────────────
extraction_qc <- bind_rows(.extraction_qc)

write_csv(pMR_all,         file.path(OUT_DIR, "pMR_all_results.csv"))
write_csv(pMR_ivw,         file.path(OUT_DIR, "pMR_ivw_results.csv"))
write_csv(pMR_sig,         file.path(OUT_DIR, "pMR_FDR_significant.csv"))
write_csv(evidence_matrix, file.path(OUT_DIR, "pMR_evidence_matrix.csv"))
write_csv(extraction_qc,   file.path(OUT_DIR, "pMR_instrument_extraction_QC.csv"))

message("\n=== Proteome-MR (local mode) complete ===")
message(sprintf("Outcomes run: %d (coeliac dropped — no local file)",
                nrow(outcomes)))
message(sprintf("FDR < 0.05 hits: %d", nrow(pMR_sig)))
message("Outputs:")
message("  pMR_all_results.csv")
message("  pMR_ivw_results.csv")
message("  pMR_FDR_significant.csv")
message("  pMR_evidence_matrix.csv")
message("  pMR_instrument_extraction_QC.csv")
message()
message("Note: Method 'IVW' = Wald ratio for proteins with only 1 instrument.")
message("      Run Script 01 (with OpenGWAS JWT) for multi-SNP IVW estimates.")
