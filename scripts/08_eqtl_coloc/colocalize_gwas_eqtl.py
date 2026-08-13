import os
import requests
import pandas as pd
import numpy as np
import scipy.stats as stats

print("=======================================================")
print("Starting formal statistical Colocalization (GWAS vs eQTL)")
print("=======================================================")

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
FILES = {
    'Endometriosis': os.path.join(OUT_DIR, 'Endometriosis_Local_Stats.csv'),
    'Fibroids': os.path.join(OUT_DIR, 'Fibroids_Local_Stats.csv')
}

GENES = ["WNT4", "GREB1", "ESR1", "WT1", "HMGA2", "PDGFRA", "STON2", "CDKN2B-AS1", "FSHB"]
REPRO_TISSUES = ["Uterus", "Ovary", "Vagina", "Cervix_Endocervix", "Cervix_Ectocervix"]

# 1. Map Genes to GENCODE
ref_url = "https://gtexportal.org/api/v2/reference/gene"
gencode_map = {}
for g in GENES:
    try:
        r = requests.get(ref_url, params={"geneId": g, "format": "json"}).json()
        if "data" in r and len(r["data"]) > 0:
            gencode_map[g] = r["data"][0]["gencodeId"]
    except: pass

# 2. Bayesian Colocalization Function (coloc-style)
def run_coloc(df1, df2, trait1_name, trait2_name):
    """
    df1: GWAS df with 'beta', 'se', 'pos'
    df2: eQTL df with 'nes', 'pValue', 'pos'
    Returns posterior probabilities for H0, H1, H2, H3, H4
    """
    # Merge on position
    merged = pd.merge(df1, df2, on='pos', suffixes=('_1', '_2'))
    if len(merged) < 10: return None
    
    # Calculate log Bayes Factors
    # For GWAS (trait 1)
    z1 = merged['beta'] / merged['se']
    v1 = merged['se']**2
    w1 = 0.2**2 # prior variance
    r1 = w1 / (v1 + w1)
    lbf1 = 0.5 * (np.log(1 - r1) + r1 * z1**2)
    
    # For eQTL (trait 2)
    # GTEx gives NES and PValue. We need SE.
    # z = abs(stats.norm.ppf(pval/2))
    # se = abs(nes / z)
    merged['z2'] = abs(stats.norm.ppf(merged['pValue'] / 2))
    merged['se2'] = abs(merged['nes'] / (merged['z2'] + 1e-9))
    z2 = merged['z2']
    v2 = merged['se2']**2
    w2 = 0.2**2 
    r2 = w2 / (v2 + w2)
    lbf2 = 0.5 * (np.log(1 - r2) + r2 * z2**2)
    
    # Hypothesis priors
    p1 = 1e-4
    p2 = 1e-4
    p12 = 1e-5
    
    # Log-sum-exp trick to handle probabilities
    def logsum(lx):
        m = np.max(lx)
        return m + np.log(np.sum(np.exp(lx - m)))
    
    # H1: Trait 1 only
    h1_lks = logsum(lbf1) + np.log(p1)
    # H2: Trait 2 only
    h2_lks = logsum(lbf2) + np.log(p2)
    # H3: Both traits, different variants
    h3_lks = logsum(lbf1) + logsum(lbf2) + np.log(p1) + np.log(p2)
    # H4: Shared variant
    h4_lks = logsum(lbf1 + lbf2) + np.log(p12)
    # H0: None
    h0_lks = 0
    
    lks = np.array([h0_lks, h1_lks, h2_lks, h3_lks, h4_lks])
    probs = np.exp(lks - logsum(lks))
    
    return {
        'PP_H0': probs[0], 'PP_H1': probs[1], 'PP_H2': probs[2], 'PP_H3': probs[3], 'PP_H4': probs[4],
        'n_snps': len(merged), 'lead_variant': merged.loc[merged['beta'].abs().idxmax()]['rsids' if 'rsids' in merged.columns else 'pos']
    }

# 3. Main Loop
coloc_results = []
# Load GWAS subset data
gwas_data = {}
for pheno, path in FILES.items():
    gwas_data[pheno] = pd.read_csv(path)
    # Convert chromosome and pos to numeric for merging
    gwas_data[pheno]['#chrom'] = pd.to_numeric(gwas_data[pheno]['#chrom'])
    gwas_data[pheno]['pos'] = pd.to_numeric(gwas_data[pheno]['pos'])

# Target Regions (from LAVA hits)
regions = [
    (1, 1376205, 2215495, "Block_2"),
    (1, 3582048, 4281476, "Block_5"),
    (1, 10753428, 11709173, "Block_13"),
    (1, 11709174, 13481025, "Block_14"),
    (1, 13481026, 14720131, "Block_15"),
    (1, 15755231, 16732168, "Block_17"),
    (1, 16732169, 17557746, "Block_18")
]

print(f"Testing {len(regions)} regions across {len(GENES)} genes and {len(REPRO_TISSUES)} tissues...")

for chrom, start, stop, block_id in regions:
    print(f"Processing {block_id}...")
    for trait_name, gwas_df in gwas_data.items():
        # Subset GWAS to this block
        gwas_subset = gwas_df[(gwas_df['#chrom'] == chrom) & (gwas_df['pos'] >= start) & (gwas_df['pos'] <= stop)].copy()
        if gwas_subset.empty: continue
        
        # Proper column names for coloc
        gwas_subset = gwas_subset.rename(columns={'beta': 'beta', 'sebeta': 'se', 'pos': 'pos'})
        
        for gene_sym, gencode in gencode_map.items():
            # Check if gene is in this block window (+/- 500kb for safety)
            # Fetch associations for this gene in GTEx for relevant tissues
            for tissue in REPRO_TISSUES:
                try:
                    url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
                    params = {"datasetId": "gtex_v8", "gencodeId": gencode, "tissueSiteDetailId": tissue, "format": "json"}
                    r = requests.get(url, params=params).json()
                    if "data" in r and len(r["data"]) > 0:
                        eqtl_df = pd.DataFrame(r["data"])
                        # Extract position from variantId (chr1_12345_A_G_b38)
                        eqtl_df['pos'] = eqtl_df['variantId'].apply(lambda x: int(x.split('_')[1]))
                        
                        # run coloc
                        res = run_coloc(gwas_subset, eqtl_df, trait_name, f"{gene_sym}_{tissue}")
                        if res:
                            res.update({
                                'Block': block_id,
                                'Trait': trait_name,
                                'Gene': gene_sym,
                                'Tissue': tissue
                            })
                            coloc_results.append(res)
                except: pass

if coloc_results:
    res_df = pd.DataFrame(coloc_results)
    res_df = res_df.sort_values('PP_H4', ascending=False)
    out_file = os.path.join(OUT_DIR, "Colocalization_GWAS_eQTL_Results.csv")
    res_df.to_csv(out_file, index=False)
    print(f"Successfully identified {len(res_df)} colocalization pairs.")
    print(f"Saved to {out_file}")
    
    # Print high-confidence hits
    hi_hits = res_df[res_df['PP_H4'] > 0.5]
    if not hi_hits.empty:
        print("\nStrong/Suggestive Colocalization Hits (PP.H4 > 0.5):")
        print(hi_hits[['Block', 'Trait', 'Gene', 'Tissue', 'PP_H4']])
else:
    print("No colocalization signals reached required snp-overlap for calculation.")

print("\nPipeline Complete!")
