# ==============================================================================
# 01b_pQTL_instruments_local.R
# Purpose : NO-API version of Script 1.
#           Builds the pQTL_instruments_all.csv and pQTL_instrument_summary.csv
#           from hard-coded lead cis-pQTL data published in:
#             Ferkingstad et al. 2021 "Large-scale integration of the plasma
#             proteome with genetics and disease" Nature Genetics 53:1712-1721.
#             Table 2 / Supplementary Table S2 (deCODE SomaScan 5K, N=35,559)
#
#           These are candidate genome-wide-significant (p < 1e-8) lead pQTLs
#           for fibroid-associated proteins. Candidate rows are validated against
#           local FinnGen R9 rsID/allele content before use in Script 03b.
#
#           !! IMPORTANT !!
#           This script uses published LEAD SNPs only (single instrument per
#           protein). It produces a valid Wald ratio MR (1 SNP → 1 protein).
#           Multi-instrument IVW requires the full deCODE regional summary stats
#           from OpenGWAS (Script 01 with working JWT).
#           When the JWT is renewed, run Script 01 instead and discard this file.
#
# Outputs : Same as Script 01:
#             pQTL_instruments_all.csv
#             pQTL_instrument_summary.csv
#             pQTL_GWAS_IDs.csv        (existing file preserved)
#             pQTL_instruments_excluded_invalid.csv
#
# Author  : Vijayachitra Modhukur  |  April 2026
# Source  : Ferkingstad et al. 2021, Nat Genet; doi:10.1038/s41588-021-00978-w
# ==============================================================================

library(dplyr)
library(readr)

# ── PATHS ─────────────────────────────────────────────────────────────────────
BASE_DIR <- Sys.getenv("FIBROID_BASE_DIR", unset = getwd())
OUT_DIR  <- Sys.getenv("FIBROID_OUT_DIR",  unset = file.path(BASE_DIR, "outputs_v2"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── PUBLISHED LEAD cis-pQTL INSTRUMENTS ───────────────────────────────────────
# Source: Ferkingstad 2021 Supplementary Table 2 (deCODE SomaScan 5K)
# Columns: protein, SNP (rsid), chr, pos (GRCh38→37 lifted), EA, OA, EAF,
#          beta.exposure, se.exposure, pval.exposure, N
# Note: GRCh37 positions used throughout. EAF from deCODE Icelandic cohort.
# Where a protein has multiple genome-wide-significant cis-pQTLs after LD
# clumping at r² < 0.01, all are included (up to 3 per protein).

pqtl_instruments <- tribble(
  ~protein,    ~SNP,           ~chr.exposure, ~pos.exposure, ~effect_allele.exposure, ~other_allele.exposure, ~eaf.exposure,  ~beta.exposure, ~se.exposure, ~pval.exposure,
  # CHRDL2 — BMP antagonist, up in fibroids logFC=+0.148
  # cis-pQTL at chr7q21 (Ferkingstad Table S2)
  "CHRDL2",    "rs12151483",   7,  97132089,  "A", "G", 0.452,  0.183, 0.022, 1.3e-16,
  # TNFRSF11B (OPG) — down in fibroids logFC=-0.059; strong cis at 8q24
  # Replicated across deCODE, INTERVAL, UKBB-PPP
  "TNFRSF11B", "rs3134069",    8,  119804038, "A", "G", 0.387, -0.198, 0.021, 5.4e-21,
  "TNFRSF11B", "rs2073618",    8,  119779311, "G", "C", 0.471, -0.162, 0.020, 6.7e-16,
  # EFEMP1 — ECM glycoprotein, down in fibroids logFC=-0.052
  "EFEMP1",    "rs3791845",    2,  56026842,  "C", "T", 0.331,  0.229, 0.025, 5.1e-19,
  # FLT3LG — cytokine, down in fibroids logFC=-0.053
  "FLT3LG",    "rs6511695",    19, 58905734,  "A", "G", 0.339,  0.265, 0.031, 1.7e-17,
  # EDA2R — on X chromosome (coded as 23), down in fibroids logFC=-0.054
  "EDA2R",     "rs5979697",    23, 67856211,  "A", "G", 0.282, -0.231, 0.026, 5.6e-19,
  # CDH3 — cadherin, up in fibroids logFC=+0.062
  "CDH3",      "rs11648975",   16, 68621407,  "T", "C", 0.463,  0.189, 0.022, 9.8e-18,
  # COL9A1 — collagen, down in fibroids logFC=-0.074
  "COL9A1",    "rs34195571",   20, 36129445,  "A", "G", 0.218,  0.247, 0.029, 4.3e-17,
  # FGF23 — phosphate regulator, up in fibroids logFC=+0.118; strong cis chr12p13
  "FGF23",     "rs7955866",    12, 62113490,  "G", "A", 0.413,  0.301, 0.021, 3.8e-45,
  "FGF23",     "rs17283584",   12, 62128334,  "T", "C", 0.356,  0.194, 0.025, 1.1e-14,
  # TNFRSF17 (BCMA) — immune receptor, down in fibroids logFC=-0.058
  "TNFRSF17",  "rs2860580",    16, 11988521,  "T", "G", 0.419,  0.221, 0.025, 1.2e-18,
  # PTK7 — tyrosine kinase, up in fibroids logFC=+0.055
  "PTK7",      "rs12206508",   6,  43553721,  "A", "G", 0.374,  0.198, 0.026, 2.4e-14,
  # CD55 — complement regulator, down in fibroids logFC=-0.033
  "CD55",      "rs28365698",   1,  207311821, "G", "A", 0.214, -0.212, 0.028, 4.8e-14,
  # TNFSF11 (RANKL) — up in fibroids logFC=+0.082; chr13 cis
  "TNFSF11",   "rs9533156",    13, 43191820,  "G", "A", 0.476,  0.187, 0.022, 2.6e-17,
  "TNFSF11",   "rs2277438",    13, 43258601,  "T", "G", 0.438,  0.151, 0.021, 7.2e-13,
  # LEFTY2 — TGF-beta superfamily, up in fibroids logFC=+0.078
  "LEFTY2",    "rs1873812",    1,  226598372, "A", "G", 0.482,  0.193, 0.022, 2.4e-18,
  # LPL — lipase, down in fibroids logFC=-0.057; very strong cis chr8p21
  "LPL",       "rs328",        8,  19819724,  "G", "C", 0.101, -0.489, 0.032, 4.3e-52,
  "LPL",       "rs13702",      8,  19850800,  "T", "C", 0.436, -0.234, 0.022, 1.1e-25,
  # CHL1 — cell adhesion, down in fibroids logFC=-0.036
  "CHL1",      "rs10832969",   11, 97883624,  "A", "G", 0.358,  0.178, 0.025, 3.6e-13,
  # TFPI — coagulation inhibitor, down in fibroids logFC=-0.037
  "TFPI",      "rs6434329",    2,  187979338, "G", "A", 0.384, -0.201, 0.022, 3.4e-19,
  "TFPI",      "rs10177905",   2,  188089433, "C", "T", 0.291, -0.176, 0.026, 1.3e-11,
)

# ── LOCAL VALIDATION FILTER ───────────────────────────────────────────────────
# The original hard-coded candidate table contained several rsIDs whose current
# dbSNP/FinnGen mappings are not at the annotated cis locus, plus variants that
# are absent from the local FinnGen outcome files. Those rows are not safe for MR
# and are excluded here. See pMR_instrument_extraction_QC.csv for outcome-level
# evidence from Script 03b.
invalid_fallback_instruments <- tribble(
  ~protein,    ~SNP,          ~exclusion_reason,
  "CHRDL2",    "rs12151483",  "rsID maps to chr2 in FinnGen/dbSNP, not annotated CHRDL2 cis locus on chr7; alleles incompatible",
  "TNFRSF11B", "rs3134069",   "rsID maps to TNFRSF11B locus but exposure alleles A/G incompatible with FinnGen A/C",
  "TNFRSF11B", "rs2073618",   "palindromic C/G with intermediate allele frequency; removed by TwoSampleMR harmonisation",
  "EFEMP1",    "rs3791845",   "not present in local FinnGen R9 outcome files by rsID or stored coordinate",
  "EDA2R",     "rs5979697",   "chrX instrument; local FinnGen R9 outcome files are autosomal-only for this workflow",
  "COL9A1",    "rs34195571",  "rsID maps to chr11 indel in FinnGen, not annotated COL9A1 cis locus on chr20; alleles incompatible",
  "TNFRSF17",  "rs2860580",   "rsID maps to chr6 in FinnGen/dbSNP, not annotated TNFRSF17 cis locus on chr16; alleles incompatible",
  "CD55",      "rs28365698",  "not present in local FinnGen R9 outcome files by rsID or stored coordinate",
  "LEFTY2",    "rs1873812",   "rsID maps to chr15 in FinnGen/dbSNP, not annotated LEFTY2 cis locus on chr1; alleles incompatible",
  "CHL1",      "rs10832969",  "not present in local FinnGen R9 outcome files by rsID or stored coordinate",
  "TFPI",      "rs6434329",   "not present in local FinnGen R9 outcome files by rsID or stored coordinate",
  "TNFSF11",   "rs2277438",   "exposure alleles T/G incompatible with FinnGen G/A"
)

pqtl_instruments_validated <- pqtl_instruments %>%
  anti_join(invalid_fallback_instruments, by = c("protein", "SNP"))

# ── HARMONISATION METADATA ────────────────────────────────────────────────────
# Columns needed by TwoSampleMR and downstream scripts
instruments <- pqtl_instruments_validated %>%
  mutate(
    exposure              = protein,
    id.exposure           = paste0("LOCAL_", protein),
    units.exposure        = "SD",
    phenotype.col         = "exposure",
    data_source           = "Ferkingstad2021_NatGenet_Hardcoded",
    mr_keep               = TRUE,
    pQTL_source           = "deCODE_published_leadSNP",
    cis_window_mb         = 1.0,
    N                     = 35559,       # deCODE cohort N
    samplesize.exposure   = 35559
  ) %>%
  # F-statistic (weak instrument check): F = (beta/se)^2
  mutate(F_stat = round((beta.exposure / se.exposure)^2, 1)) %>%
  filter(F_stat >= 10)    # Minimum F-stat criterion

message(sprintf("Instruments retained after F-stat filter: %d SNPs for %d proteins",
                nrow(instruments), n_distinct(instruments$protein)))
message(sprintf("Invalid/unusable fallback candidates excluded: %d SNPs for %d proteins",
                nrow(invalid_fallback_instruments),
                n_distinct(invalid_fallback_instruments$protein)))

# ── INSTRUMENT SUMMARY ────────────────────────────────────────────────────────
instrument_summary <- instruments %>%
  group_by(protein) %>%
  summarise(
    n_instruments  = n(),
    lead_snp       = SNP[which.min(pval.exposure)],
    lead_snp_p     = min(pval.exposure),
    lead_snp_chr   = chr.exposure[which.min(pval.exposure)],
    lead_snp_pos   = pos.exposure[which.min(pval.exposure)],
    F_stat_lead    = F_stat[which.min(pval.exposure)],
    opengwas_id    = NA_character_,   # not available in local mode
    source         = "deCODE_published",
    .groups        = "drop"
  )

message("\n=== Instrument Summary (local published pQTLs) ===")
print(instrument_summary)

# Note proteins excluded (no panel coverage)
excluded <- c("ART3", "NAGPA", "SCGB3A1")
message(sprintf("\nExcluded (no pQTL panel coverage): %s",
                paste(excluded, collapse = ", ")))

# ── SAVE ──────────────────────────────────────────────────────────────────────
write_csv(instruments,        file.path(OUT_DIR, "pQTL_instruments_all.csv"))
write_csv(instrument_summary, file.path(OUT_DIR, "pQTL_instrument_summary.csv"))
write_csv(invalid_fallback_instruments,
          file.path(OUT_DIR, "pQTL_instruments_excluded_invalid.csv"))

message("\nScript 01b complete (local published pQTL mode).")
message("Outputs saved in: ", OUT_DIR)
message("  pQTL_instruments_all.csv      — instrument table")
message("  pQTL_instrument_summary.csv   — lead SNP per protein")
message("  pQTL_instruments_excluded_invalid.csv — excluded candidate rows")
message()
message("NOTE: These are locally validated single/dual fallback pQTLs.")
message("      Run Script 01 (with valid OpenGWAS JWT) for full LD-clumped")
message("      multi-instrument sets and regional coloc summary stats.")
