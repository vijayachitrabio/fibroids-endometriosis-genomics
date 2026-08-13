# ==============================================================================
# run_temporal_analysis.R
# Final 'Strong Paper' Version (v2.4) - Optimized & Reviewer-Proof
# Fast array/string alignment for UKB HES variables + Advanced Metrics
# ==============================================================================

library(data.table)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)

OUT_DIR <- "outputs_0404"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

cat("\n=== 1. Fast Extraction & Alignment of Hospital Events ===\n")
dt <- fread("uf_data_0404.csv", nThread=8)

# Function to extract and align a code array to a date array
extract_pair <- function(dt, code_pref, date_pref) {
    code_cols <- grep(paste0("^", code_pref), names(dt), value=TRUE)
    date_cols <- grep(paste0("^", date_pref), names(dt), value=TRUE)
    
    if(length(code_cols) == 0 || length(date_cols) == 0) return(NULL)
    
    cat(sprintf("  Pairing %s (%d cols) -> %s (%d cols)...\n", code_pref, length(code_cols), date_pref, length(date_cols)))
    
    # Process Dates (always array)
    dt_d <- melt(dt[, c("eid", date_cols), with=FALSE], id.vars="eid", measure.vars=date_cols, value.name="date", na.rm=TRUE)
    dt_d[, idx := as.integer(str_extract(variable, "[0-9]+$"))]
    dt_d[is.na(idx), idx := 0]
    setorder(dt_d, eid, idx)
    dt_d[, rn := seq_len(.N), by=eid]

    # Process Codes
    if(length(code_cols) == 1) {
        # Concatenated string
        dt_c <- dt[, c("eid", "p21003_i3", code_cols), with=FALSE]
        setnames(dt_c, code_cols, "code_str")
        dt_c <- dt_c[code_str != ""]
        dt_c_long <- dt_c[, .(code = unlist(strsplit(as.character(code_str), "\\|"))), by=.(eid, p21003_i3)]
        dt_c_long[, rn := seq_len(.N), by=eid]
        res <- dt_c_long[dt_d, on=.(eid, rn), nomatch=0]
    } else {
        # Array columns
        dt_c <- melt(dt[, c("eid", "p21003_i3", code_cols), with=FALSE], id.vars=c("eid", "p21003_i3"), measure.vars=code_cols, value.name="code", na.rm=TRUE)
        dt_c <- dt_c[code != ""]
        dt_c[, idx := as.integer(str_extract(variable, "[0-9]+$"))]
        dt_c[is.na(idx), idx := 0]
        setorder(dt_c, eid, idx)
        res <- dt_c[dt_d, on=.(eid, idx), nomatch=0]
    }
    
    return(res[, .(eid, p21003_i3, code, date)])
}

# The defined mappings for UKB summary & split arrays
events <- list()
events[["icd10"]]  <- extract_pair(dt, "p41270", "p41280")
events[["icd9_m"]] <- extract_pair(dt, "p41202", "p41203")
events[["icd9_s"]] <- extract_pair(dt, "p41204", "p41205")

# Safely bind available results
dt_events <- rbindlist(Filter(function(x) !is.null(x) && nrow(x) > 0, events), fill=TRUE)
dt_events[, date := as.IDate(date)]
dt_events[, code := str_remove_all(as.character(code), "[^A-Z0-9]")]
cat("Total correctly paired events:", nrow(dt_events), "\n")

# --- 2. Define Traits ---
cat("\n=== 2. Trait Identification ===\n")
traits <- list(
  fibroids = "^D25", endometriosis = "^N80", heavy_menstrual_bleeding = "^N92",
  chronic_pelvic_pain = "^R102", pcos = "^E282", endometrial_polyp = "^N840",
  ibs = "^K58", fibromyalgia = "^M797", anxiety = "^F41", depression = "^F32|^F33",
  ptsd = "^F431", nafld = "^K760", vitamin_d_deficiency = "^E559",
  osteoporosis = "^M81", osteoarthritis = "^M15|^M16|^M17|^M18|^M19",
  iron_deficiency = "^D50", adenomyosis = "^N800"
)

dt_events[, trait := NA_character_]
for (tName in names(traits)) {
  dt_events[grepl(traits[[tName]], code), trait := tName]
}

# Get earliest date for each trait per person
early_events <- dt_events[!is.na(trait) & !is.na(date), .(date = min(date)), by = .(eid, trait)]

# --- 3. Advanced Temporal Summaries ---
cat("\n=== 3. Temporal Architecture (Reviewer-Proof Metrics) ===\n")

summarize_refined <- function(all_events, primary_disease, pheno_label) {
  disease_dates <- all_events[trait == primary_disease, .(eid, prim_date = date)]
  comorb_dates  <- all_events[trait != primary_disease]
  if (nrow(disease_dates) == 0) return(NULL)
  
  merged <- merge(comorb_dates, disease_dates, by="eid")
  merged[, delta_days := as.numeric(date - prim_date)]
  merged[, timing := fcase(delta_days < 0, "Before", delta_days > 0, "After", delta_days == 0, "Same Day")]
  
  sum_tab <- merged[, .(
    N_total = .N, 
    N_before = sum(timing == "Before"), 
    N_after = sum(timing == "After"),
    N_same = sum(timing == "Same Day"), 
    N_inf = sum(timing %in% c("Before", "After")),
    Med_Lag_Post_Yrs = round(median(delta_days[delta_days > 0]/365.25, na.rm=T), 1),
    Med_Lag_Pre_Yrs = round(median(abs(delta_days[delta_days < 0])/365.25, na.rm=T), 1)
  ), by = trait]
  
  sum_tab[, Prop_After_Conservative := N_after / N_total]
  
  # Flags & Decisions
  sum_tab[, Flag_High_SameDay := (N_same / N_total) > 0.5]
  sum_tab[, Flag_Low_Inf      := N_inf < 10]
  sum_tab[, Direction := fcase(
      Prop_After_Conservative > 0.6, "Predominantly Post-Disease", 
      Prop_After_Conservative < 0.4, "Predominantly Pre-Disease", 
      default = "Mixed/Bidirectional")]
  sum_tab[, Disease := pheno_label]
  return(sum_tab)
}

res_fib  <- summarize_refined(early_events, "fibroids", "Fibroids")
res_endo <- summarize_refined(early_events, "endometriosis", "Endometriosis")
all_results <- rbind(res_fib, res_endo)

# --- 4. Clinical Context (Age at Diagnosis) ---
get_diag_age <- function(dt_e, pattern) {
    ages <- dt_e[grepl(pattern, code), p21003_i3]
    return(median(ages, na.rm=T))
}
age_fib <- get_diag_age(dt_events, "^D25")
age_endo <- get_diag_age(dt_events, "^N80")

# --- 5. Save Results ---
if (!is.null(all_results)) {
    fwrite(all_results, file.path(OUT_DIR, "Table_Temporal_Advanced_Metrics.csv"))
    
    # Numeric text summary for Manuscript (Apply N >= 20 filter)
    filtered <- all_results[N_total >= 20]
    text_summary <- filtered[, .(
        Traits = .N, 
        Post = sum(Prop_After_Conservative > 0.6), 
        Pre = sum(Prop_After_Conservative < 0.4), 
        Med_Prop = round(median(Prop_After_Conservative, na.rm=T),3)
    ), by=Disease]
    
    fwrite(text_summary, file.path(OUT_DIR, "Table_Temporal_Text_Summary.csv"))
    print(text_summary)
    cat(sprintf("\n[Clinical Context]\n- Fibroids Med Age: %.1f\n- Endometriosis Med Age: %.1f\n", age_fib, age_endo))

    # --- 6. Plotting ---
    plot_data <- filtered %>% mutate(Trait_Label = str_to_title(str_replace_all(trait, "_", " ")))
    p <- ggplot(plot_data, aes(x = Prop_After_Conservative, y = reorder(Trait_Label, Prop_After_Conservative))) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
      geom_segment(aes(x = 0, xend = Prop_After_Conservative, yend = Trait_Label), color = "grey85") +
      geom_point(aes(color = Direction, size = N_total, shape = Flag_High_SameDay)) +
      facet_wrap(~Disease) +
      theme_minimal() +
      scale_color_manual(values = c("Predominantly Post-Disease" = "#2C3E50", "Predominantly Pre-Disease" = "#E74C3C", "Mixed/Bidirectional" = "#95A5A6")) +
      scale_shape_manual(name = "High Same-Day (>50%)", values = c("FALSE" = 19, "TRUE" = 4)) +
      labs(title = "Temporal Architecture: Disease vs. Comorbidity", 
           subtitle = "Filter N>=20. Conservative Metric. High same-day denoted by crosses.",
           x = "Proportion AFTER Diagnosis", y = NULL) +
      theme(legend.position = "bottom", axis.text.y = element_text(size = 8), legend.direction = "vertical")
    
    ggsave(file.path(OUT_DIR, "Figure_Temporal_Reviewer_Proof_Final.png"), p, width = 11, height = 9, dpi = 300)
    cat("\nAnalysis complete! Results saved in", OUT_DIR, "\n")
} else {
    cat("\nERROR: No events recovered.\n")
}
