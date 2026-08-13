import sys, os, gzip, subprocess
import pandas as pd
import numpy as np
from scipy import stats

# Configuration
import os
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
GWAS_DIR = os.path.join(BASE_DIR, "Finngen_phewas", "finngen_r9_phewas")
OUTPUT_RES = "MR_PheWAS_Robust_Results.csv"

# Full 44-trait catalogue (R9)
FINNGEN_CATALOGUE = {
    'irritable_bowel_syndrome':   ('K11_IBS',                       True),
    'endometrial_cancer':         ('C3_CORPUS_UTERI_EXALLC',        True),
    'depression':                 ('F5_DEPRESSIO',                   False),
    'anxiety':                    ('F5_ANXIETY',                     False),
    'heavy_menstrual_bleeding':   ('N14_MENORRHAGIA',                False),
    'pcos':                       ('E4_PCOS',                        False),
    'ovarian_cyst':               ('N14_OVARYCYST',                  False),
    'ovarian_cancer':             ('C3_OVARY_EXALLC',                False),
    'breast_cancer':              ('C3_BREAST_EXALLC',                False),
    'cervical_cancer':            ('C3_CERVIX_UTERI_EXALLC',         False),
    'endometrial_polyp':          ('N14_UTER_POLYP',                 False),
    'adenomyosis':                ('N14_ADENOMYOSIS_UTERUS',         False),
    'obesity':                    ('E4_OBESITY',                     False),
    'type2_diabetes':             ('T2D',                            False),
    'hypertension':               ('I9_HYPTENS',                     False),
    'hypothyroidism':             ('E4_HYTHY_AUTO',                  False),
    'hyperthyroidism':            ('AUTOIMMUNE_HYPERTHYROIDISM',     False),
    'osteoarthritis':             ('M13_ARTHROSIS',                  False),
    'osteoporosis':               ('M13_OSTEOPOROSIS',               False),
    'rheumatoid_arthritis':       ('M13_RHEUMA',                     False),
    'sle':                        ('M13_SYSTSLUPUSERYT',             False),
    'ibd':                        ('K11_IBD_STRICT',                 False),
    'coeliac_disease':            ('K11_COELIAC',                    False),
    'migraine':                   ('G6_MIGRAINE',                    False),
    'fibromyalgia':               ('M13_FIBROMYALGIA',               False),
    'chronic_pelvic_pain':        ('N14_FEMGENPAIN',                 False),
    'iron_deficiency_anaemia':    ('D3_ANAEMIA_IRONDEF',             False),
    'recurrent_uti':              ('N14_URINARYINFECT_RECUR',        False),
    'overactive_bladder':         ('N14_NEUROMUSCDYSBLADD',          False),
    'urinary_incontinence':       ('N14_STRESSINCONT',               False),
    'dyslipidaemia':              ('E4_LIPOPROT',                    False),
    'nafld':                      ('NAFLD',                          False),
    'gerd':                       ('K11_REFLUX',                     False),
    'peptic_ulcer_disease':       ('K11_PEPTICULCER',                False),
    'back_pain':                  ('M13_LOWBACKPAIN',                False),
    'psoriasis':                  ('L12_PSORIASIS',                  False),
    'vitamin_d_deficiency':       ('E4_VIT_D_DEF',                   False),
    'leiomyosarcoma':             ('C3_LEIOMYOSARCOMA_EXALLC',       False),
    'interstitial_cystitis':      ('INTERSTITIAL_CYSTIT_CHRONIC',    False),
    'vulvodynia':                 ('N14_VULVODYN',                   False),
    'thrombocytopenia':           ('D3_THROMBOCYTOPENIANAS',         False),
    'ptsd':                       ('F5_PTSD',                        False),
    'endometriosis':              ('N14_ENDOMETRIOSIS',              True),
    'uterine_fibroids':           ('CD2_BENIGN_LEIOMYOMA_UTERI',     True),
}

# Exposure-specific P-value thresholds (relax for smaller GWAS)
EXPOSURE_THRESHOLDS = {
    'Endometriosis':     5e-8,
    'Uterine_Fibroids':  5e-8,
    'Endometrial_CA':    1e-5,
}

EXPOSURE_FILES = {
    'Endometriosis':     os.path.join(GWAS_DIR, 'finngen_R9_N14_ENDOMETRIOSIS.gz'),
    'Uterine_Fibroids':  os.path.join(GWAS_DIR, 'finngen_R9_CD2_BENIGN_LEIOMYOMA_UTERI.gz'),
    'Endometrial_CA':    os.path.join(GWAS_DIR, 'finngen_R9_C3_CORPUS_UTERI_EXALLC.gz'),
}

def get_available_outcomes():
    available = {}
    files = os.listdir(GWAS_DIR)
    for name, (code, is_expo) in FINNGEN_CATALOGUE.items():
        # Check for prefixed or raw phenocode
        found = None
        if f"finngen_R9_{code}.gz" in files: found = f"finngen_R9_{code}.gz"
        elif f"{code}.gz" in files: found = f"{code}.gz"
        
        if found:
            path = os.path.join(GWAS_DIR, found)
            if os.path.getsize(path) > 1024*1024: # > 1MB
                available[name] = path
    return available

def read_sumstats(file_path, name, is_exposure=True):
    print(f"  Reading {os.path.basename(file_path)} ...")
    try:
        df = pd.read_csv(file_path, sep='\t', compression='gzip')
    except Exception as e:
        print(f"    ✗ Error reading {name}: {e}")
        return pd.DataFrame()

    # Column Mapping
    col_map = {'#chrom':'chr', 'pos':'pos', 'ref':'oa', 'alt':'ea', 'rsids':'rsid', 'pval':'pval', 'beta':'beta', 'sebeta':'se'}
    df = df.rename(columns=col_map)
    df = df[['chr','pos','ea','oa','rsid','pval','beta','se']]

    if is_exposure:
        # P-value Threshold
        threshold = EXPOSURE_THRESHOLDS.get(name, 5e-8)
        exp_hits = df[df['pval'] < threshold].copy()
        if exp_hits.empty:
            print(f"  ✗ No instruments found at P < {threshold}")
            return pd.DataFrame()
        
        # Simple LD Pruning (500kb)
        exp_hits = exp_hits.sort_values('pval')
        pruned = []
        for _, row in exp_hits.iterrows():
            if not any((abs(row['pos'] - p['pos']) < 500000 and row['chr'] == p['chr']) for p in pruned):
                pruned.append(row)
        
        print(f"    → {len(pruned)} instruments (P<{threshold}, LD-pruned 500kb)")
        return pd.DataFrame(pruned)
    return df

def rc(allele):
    mapping = str.maketrans("ATCG", "TAGC")
    return allele.upper().translate(mapping)

def lookup_outcome(instruments, outcome_file):
    print(f"  Looking up in {os.path.basename(outcome_file)} ...")
    snp_set = set(instruments['rsid'].dropna().tolist())
    rows = []
    with gzip.open(outcome_file, 'rt') as f:
        header_line = f.readline().lstrip('#').strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if not parts or len(parts) < 7: continue
            rsid_field = parts[4]
            for rsid in rsid_field.split(','):
                if rsid.strip() in snp_set:
                    rows.append(dict(zip(header_line, parts)))
                    break
    
    if not rows:
        print("    → 0 SNPs matched outcome")
        return pd.DataFrame()
        
    df = pd.DataFrame(rows)
    df.rename(columns={'sebeta':'se_out','rsids':'rsid','chrom':'chr',
                        'af_alt':'eaf_out','ref':'oa_out','alt':'ea_out', 'beta':'beta_out', 'pval':'pval_out'}, inplace=True)
    df['beta_out'] = pd.to_numeric(df['beta_out'], errors='coerce')
    df['se_out'] = pd.to_numeric(df['se_out'], errors='coerce')
    df['pval_out'] = pd.to_numeric(df['pval_out'], errors='coerce')
    df = df[['rsid','ea_out','oa_out','beta_out','se_out','pval_out']].copy()
    print(f"    → {len(df)} SNPs matched in outcome")
    return df

def mr_ivw(b_exp, se_exp, b_out, se_out):
    if len(b_exp) < 2: return np.nan, np.nan, np.nan
    w = 1 / (se_out**2)
    beta_ivw = np.sum(w * b_exp * b_out) / np.sum(w * b_exp**2)
    se_ivw = np.sqrt(1 / np.sum(w * b_exp**2))
    pval = 2 * (1 - stats.norm.cdf(abs(beta_ivw / se_ivw)))
    return beta_ivw, se_ivw, pval

def mr_weighted_median(b_exp, se_exp, b_out, se_out):
    if len(b_exp) < 3: return np.nan, np.nan, np.nan
    # Convert to numpy to avoid index mismatches
    b_exp_v, se_exp_v = b_exp.values, se_exp.values
    b_out_v, se_out_v = b_out.values, se_out.values
    
    b_ratio = b_out_v / b_exp_v
    w = (b_exp_v**2) / (se_out_v**2)
    idx = np.argsort(b_ratio)
    b_ratio, w = b_ratio[idx], w[idx]
    cum_w = np.cumsum(w) / np.sum(w)
    median_beta = np.interp(0.5, cum_w, b_ratio)
    # Simple SE via bootstrapping (simplified)
    se = np.nan 
    pval = np.nan 
    return median_beta, se, pval

def steiger_filter(df):
    # Convert to numpy for calculations to avoid index issues
    b_exp, se_exp = df['beta'].values, df['se'].values
    b_out, se_out = df['beta_out'].values, df['se_out'].values
    
    z_exp = abs(b_exp / se_exp)
    z_out = abs(b_out / se_out)
    ok = z_exp > (z_out + 0.1)
    
    n_filtered = len(df) - sum(ok)
    if n_filtered > 0:
        print(f"    → Steiger filtered {n_filtered} SNPs (suggesting reverse causality)")
    return df.iloc[np.where(ok)[0]].copy()

def run_phewas_all():
    print("="*60)
    print("PheWAS-MR Analysis (Robustly Corrected) — FinnGen R9")
    print("="*60)
    
    outcomes = get_available_outcomes()
    results = []
    
    for exp_name, exp_path in EXPOSURE_FILES.items():
        print(f"\n▶ Exposure: {exp_name}")
        exp_instruments = read_sumstats(exp_path, exp_name, is_exposure=True)
        if exp_instruments.empty: continue
        
        for out_name, out_path in outcomes.items():
            # Testing everything including core exposures
            print(f"  ▶ Outcome: {out_name}")
            out_data = lookup_outcome(exp_instruments, out_path)
            if out_data.empty: continue
            
            # Harmonise
            merged = exp_instruments.merge(out_data, on='rsid', how='inner')
            if merged.empty: continue
            
            # Steiger Filter
            merged = steiger_filter(merged)
            if len(merged) < 1:
                print("    ✗ Only 0 harmonised SNPs remain after filtering — skipping")
                continue
                
            # Wald Ratio for 1 SNP
            if len(merged) == 1:
                beta = merged['beta_out'].values[0] / merged['beta'].values[0]
                se = merged['se_out'].values[0] / abs(merged['beta'].values[0])
                pval = merged['pval_out'].values[0]
                res = {'Exposure': exp_name, 'Outcome': out_name, 'N_instruments': 1,
                       'IVW_OR': np.exp(beta), 'IVW_P': pval, 'WM_OR': np.nan, 'Conclusion': 'Wald Ratio'}
            else:
                beta, se, pval = mr_ivw(merged['beta'], merged['se'], merged['beta_out'], merged['se_out'])
                wm_beta, _, _ = mr_weighted_median(merged['beta'], merged['se'], merged['beta_out'], merged['se_out'])
                
                conclusion = "Null"
                if pval < 0.05: conclusion = "Nominal"
                if pval < 0.05 / len(outcomes): conclusion = "High Confidence"
                
                res = {
                    'Exposure': exp_name, 'Outcome': out_name, 'N_instruments': len(merged),
                    'IVW_OR': np.exp(beta), 'IVW_P': pval, 'WM_OR': np.exp(wm_beta),
                    'Conclusion': conclusion
                }
            
            results.append(res)
            print(f"    IVW OR={res['IVW_OR']:.3f}, P={res['IVW_P']:.4f} | Conclusion: {conclusion}")

    # FDR Correction
    final_df = pd.DataFrame(results)
    if not final_df.empty:
        final_df['IVW_FDR_P'] = stats.false_discovery_control(final_df['IVW_P'])
        final_df.to_csv(OUTPUT_RES, index=False)
        print(f"\nFinal Results saved to {OUTPUT_RES}")

if __name__ == "__main__":
    run_phewas_all()
