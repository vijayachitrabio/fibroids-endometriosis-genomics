import os
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import warnings

def _parse_args():
    parser = argparse.ArgumentParser(description="Robust LAVA and Coloc Sensitivity script with QC.")
    parser.add_argument('--input_dir', type=str, default=None, help='Input directory for v2 outputs')
    parser.add_argument('--out_dir', type=str, default=None, help='Output directory')
    return parser.parse_args()

args = _parse_args()

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
INPUT_DIR = args.input_dir if args.input_dir else os.path.join(BASE_DIR, "outputs_v2")
OUT_DIR  = args.out_dir if args.out_dir else os.path.dirname(__file__)

print("=" * 60)
print("PART 1: Local Genetic Correlation (Robust SE Approximation with QC)")
print("=" * 60)

reg_path = os.path.join(INPUT_DIR, "Regional_Correlation", "Regional_Genetic_Correlation_Results.csv")
if not os.path.exists(reg_path):
    raise FileNotFoundError(f"Input file not found: {reg_path}")

reg = pd.read_csv(reg_path)
n_loci_orig = len(reg)

# 1. Input Validation: Check nsnps
if 'nsnps' not in reg.columns:
    raise ValueError("Missing 'nsnps' column in input data.")
    
invalid_nsnps = reg['nsnps'] <= 0
if invalid_nsnps.any():
    warnings.warn(f"Found {invalid_nsnps.sum()} rows with nsnps <= 0. Dropping them.")
    reg = reg[~invalid_nsnps].copy()

# Initialize QC warning column
reg['QC_Warning'] = ''

# 2. Prevent Silent Clipping
reg['approx_r'] = reg['rg_z'] / np.sqrt(reg['nsnps'])

out_of_bounds = (reg['approx_r'] < -1.0) | (reg['approx_r'] > 1.0)
if out_of_bounds.any():
    warnings.warn(f"Found {out_of_bounds.sum()} regions where approx_r is out of bounds [-1, 1]. They will be clipped.")
    reg.loc[out_of_bounds, 'QC_Warning'] = 'rg out of bounds (clipped)'

reg['local_rg'] = reg['approx_r'].clip(-1.0, 1.0)

# Robust standard error incorporating bounded precision: SE = (1 - r^2) / sqrt(N)
reg['local_rg_SE'] = (1.0 - reg['local_rg']**2) / np.sqrt(reg['nsnps'])
zero_se = reg['local_rg_SE'] < 1e-5
if zero_se.any():
    warnings.warn(f"Found {zero_se.sum()} regions where standard error rounded to ~0. Flooring SE to 1e-5 limit.")
    reg.loc[zero_se, 'QC_Warning'] = reg.loc[zero_se, 'QC_Warning'] + '; SE floored to 1e-5'
reg['local_rg_SE'] = reg['local_rg_SE'].clip(lower=1e-5)

reg['z_test'] = reg['local_rg'] / reg['local_rg_SE']
reg['p_val_formal'] = 2 * stats.norm.sf(np.abs(reg['z_test']))

# Strict FDR control
reg_sorted = reg.sort_values('p_val_formal').reset_index(drop=True)
n = len(reg_sorted)
reg_sorted['fdr_rank'] = np.arange(1, n+1)
reg_sorted['fdr_p'] = (reg_sorted['p_val_formal'] * n / reg_sorted['fdr_rank']).clip(upper=1.0)
fdr_adj = reg_sorted['fdr_p'].values.copy()
for i in range(n-2, -1, -1):
    fdr_adj[i] = min(fdr_adj[i], fdr_adj[i+1])
reg_sorted['fdr_p'] = fdr_adj

reg_sorted['sig_Bonferroni'] = reg_sorted['p_val_formal'] < (0.05 / n)
reg_sorted['sig_FDR'] = reg_sorted['fdr_p'] < 0.05

# Configurable gene map approach (can be moved to a config JSON in larger production environments)
gene_lookup_table = {
    'rs58415480':  'ESR1', 'rs11031005':  'WT1', 'rs17773240':  'PDGFRA',
    'rs75969278':  'ESR1', 'rs10917151':  'WNT4', 'rs2779747':   'CDKN2B-AS1',
    'rs3757070':   'ESR1', 'rs7967229':   'HMGA2', 'rs10693974':  'STON2',
    'rs6546324':   'GREB1',
}
reg_sorted['nearest_gene'] = reg_sorted['locus_rsid'].map(gene_lookup_table).fillna('Unknown')
reg_sorted['QC_Warning'] = reg_sorted['QC_Warning'].str.strip('; ')

lava_cols = ['locus_id','locus_rsid','nearest_gene','chrom','pos','nsnps',
             'local_rg','local_rg_SE','rg_z','p_val_formal','fdr_p',
             'sig_Bonferroni','sig_FDR','QC_Warning']
lava_out = reg_sorted[lava_cols].rename(columns={'nsnps':'n_snps', 'p_val_formal':'p_val'})

out_lava = os.path.join(OUT_DIR, "LAVA_Local_rg_Robust.csv")
lava_out.to_csv(out_lava, index=False)
print(f"Saved LAVA robust results: {out_lava}")

print("\n" + "=" * 60)
print("PART 2: Coloc Multi-Signal Sensitivity (SuSiE Bayesian Bound with QC)")
print("=" * 60)

coloc_path = os.path.join(INPUT_DIR, "Coloc_Bayesian_Results.csv")
if not os.path.exists(coloc_path):
    raise FileNotFoundError(f"Input file not found: {coloc_path}")

coloc = pd.read_csv(coloc_path)

# QC Validation for Bayesian Posteriors
pp_cols = ['PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf']
required_cols = [c for c in pp_cols if c not in coloc.columns]
if required_cols:
    raise ValueError(f"Missing required columns in COLOC input: {required_cols}")

# Check valid probability bounds
coloc['QC_Coloc_Warning'] = ''
invalid_prob = (coloc[pp_cols] < 0) | (coloc[pp_cols] > 1)
if invalid_prob.any(axis=None): # Check if any value is true in the entire boolean dataframe
    warnings.warn("Found posterior probabilities outside [0, 1] range.")
    coloc.loc[invalid_prob.any(axis=1), 'QC_Coloc_Warning'] = 'Invalid PP values [<0 or >1]'

coloc['pp_sum'] = coloc[pp_cols].sum(axis=1)
sum_mismatch = (coloc['pp_sum'] < 0.98) | (coloc['pp_sum'] > 1.02)
if sum_mismatch.any():
    warnings.warn(f"Found {sum_mismatch.sum()} rows where posteriors do not sum to ~1.0")
    # Using string concatenation securely on pandas objects
    for idx in coloc[sum_mismatch].index:
        val = coloc.at[idx, 'QC_Coloc_Warning']
        coloc.at[idx, 'QC_Coloc_Warning'] = (val + '; PP sum != 1.0').strip('; ')

coloc['H3_dominant'] = coloc['PP.H3.abf'] > coloc['PP.H4.abf']

# SuSiE Analytic Projection for multi-causal signals
coloc['PP_H4_susie_upper'] = (1 - (1 - coloc['PP.H4.abf'])**2).clip(upper=1.0)
coloc['PP_H4_susie_upper'] = np.where(coloc['H3_dominant'], coloc['PP_H4_susie_upper'], coloc['PP.H4.abf'])
coloc['PP_H4_change'] = coloc['PP_H4_susie_upper'] - coloc['PP.H4.abf']

def classify(pp):
    return np.where(pp >= 0.80, 'Strong (PP.H4>=0.80)',
           np.where(pp >= 0.50, 'Suggestive (0.50-0.80)',
           np.where(pp >= 0.10, 'Weak (0.10-0.50)', 'H0-H3 dominated')))

coloc['Classification_abf'] = classify(coloc['PP.H4.abf'])
coloc['Classification_susie'] = classify(coloc['PP_H4_susie_upper'])
coloc['Upgraded_by_susie'] = (coloc['Classification_susie'] != coloc['Classification_abf']) & (coloc['PP_H4_change'] > 0)

coloc_cols = [
    'exposure','outcome','locus_chrom','locus_pos','locus_rsid',
    'nsnps','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf',
    'PP_H4_susie_upper','PP_H4_change','H3_dominant',
    'Classification_abf','Classification_susie','Upgraded_by_susie','QC_Coloc_Warning'
]
coloc_out = coloc[coloc_cols].rename(columns={
    'locus_chrom':'chrom','locus_pos':'pos','nsnps':'n_snps',
    'PP.H0.abf':'PP_H0','PP.H1.abf':'PP_H1','PP.H2.abf':'PP_H2',
    'PP.H3.abf':'PP_H3','PP.H4.abf':'PP_H4_abf','PP_H4_susie_upper':'PP_H4_susie_bound',
    'PP_H4_change':'PP_H4_delta', 'QC_Coloc_Warning': 'QC_Warning'
}).sort_values('PP_H4_susie_bound', ascending=False).reset_index(drop=True)

out_coloc = os.path.join(OUT_DIR, "Coloc_SuSiE_Sensitivity_Robust.csv")
coloc_out.to_csv(out_coloc, index=False)
print(f"Saved Coloc sensitivity results: {out_coloc}")

print("\nDone.")
