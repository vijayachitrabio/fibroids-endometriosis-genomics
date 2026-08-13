import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import os

def create_publication_plots(data_path, output_dir):
    # Set modern publication-ready theme
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.linewidth": 1.2,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight"
    })
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Load Data
    print(f"Loading data from {data_path}...")
    df = pd.read_csv(data_path)
    
    # Filter out valid rows
    df = df.dropna(subset=['local_rg', 'p_val'])
    
    # Calculate -log10 P-value and Confidence Intervals
    df['-log10(P)'] = -np.log10(df['p_val'])
    df['rg_lower'] = df['local_rg'] - 1.96 * df['local_rg_SE']
    df['rg_upper'] = df['local_rg'] + 1.96 * df['local_rg_SE']
    
    # Sort by significance
    df = df.sort_values('p_val').reset_index(drop=True)
    
    # Get top 20 for Forest Plot
    top20 = df.head(20).copy()
    
    # --- PLOT 1: FOREST PLOT ---
    print("Generating Forest Plot...")
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Colors for positive/negative
    colors = ['#d73027' if rg > 0 else '#4575b4' for rg in top20['local_rg']]
    
    # Define Y-axis positions (reverse so rank 1 is at top)
    y_pos = np.arange(len(top20))[::-1]
    
    # Plot Confidence Intervals
    ax.hlines(y_pos, top20['rg_lower'], top20['rg_upper'], color=colors, lw=2, alpha=0.7)
    
    # Plot Point Estimates
    sizes = np.clip(100 * (1 / (top20['local_rg_SE'] + 0.01)), 50, 300)
    ax.scatter(top20['local_rg'], y_pos, color=colors, s=sizes, edgecolor='white', lw=0.5, zorder=3)
    
    # Vertical line at 0
    ax.axvline(0, color='black', lw=1, linestyle='--', alpha=0.5, zorder=1)
    
    # Y-axis Labels (Gene or Locus ID)
    labels = top20['nearest_gene'].fillna('Unknown').astype(str)
    labels = [f"{gene} (Chr{chrom})" if gene != 'Unknown' else f"Locus Chr{chrom}" for gene, chrom in zip(labels, top20['chrom'])]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontweight='bold')
    
    # Titles and Labels
    ax.set_xlabel('Local Genetic Correlation ($r_g$)', fontweight='bold')
    ax.set_title('Top 20 Local Genetic Correlations\nEndometriosis vs. Uterine Fibroids', fontweight='bold', pad=15)
    
    # Strip unnecessary spines
    sns.despine(left=True)
    ax.xaxis.grid(True, linestyle=':', alpha=0.6)
    
    plt.savefig(os.path.join(output_dir, 'Forest_Plot_Top_Loci.pdf'))
    plt.savefig(os.path.join(output_dir, 'Forest_Plot_Top_Loci.png'))
    plt.close()
    
    # --- PLOT 2: VOLCANO PLOT ---
    print("Generating Volcano Plot...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sig_thresh = -np.log10(0.05 / 2495)  # Bonferroni approx
    
    # Default color
    point_colors = np.where((df['-log10(P)'] > sig_thresh) & (df['local_rg'] > 0), '#d73027',
                    np.where((df['-log10(P)'] > sig_thresh) & (df['local_rg'] < 0), '#4575b4', '#bdbdbd'))
    point_alphas = np.where(df['-log10(P)'] > sig_thresh, 0.8, 0.3)
    point_sizes = np.where(df['-log10(P)'] > sig_thresh, 60, 20)
    
    ax.scatter(df['local_rg'], df['-log10(P)'], c=point_colors, s=point_sizes, alpha=0.7, edgecolor='white', lw=0.3)
    
    ax.axhline(sig_thresh, color='black', lw=1, linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', lw=1, linestyle='-', alpha=0.3)
    
    # Annotate top hits
    for i, row in top20.head(10).iterrows():
        label = str(row['nearest_gene']) if str(row['nearest_gene']) != 'Unknown' else f"Chr{row['chrom']}"
        ax.annotate(label, (row['local_rg'], row['-log10(P)']),
                   xytext=(5, 5), textcoords='offset points', fontsize=8,
                   bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.7))
    
    ax.set_xlabel('Local Genetic Correlation ($r_g$)', fontweight='bold')
    ax.set_ylabel('-log$_{10}$ (p-value)', fontweight='bold')
    ax.set_title('Local Genetic Correlation Volcano Plot', fontweight='bold', pad=15)
    sns.despine()
    
    plt.savefig(os.path.join(output_dir, 'Volcano_Plot_LAVA.pdf'))
    plt.savefig(os.path.join(output_dir, 'Volcano_Plot_LAVA.png'))
    plt.close()
    
    # --- PLOT 3: LOCAL MANHATTAN PLOT ---
    print("Generating Manhattan Plot...")
    fig, ax = plt.subplots(figsize=(12, 5))
    
    colors_manh = ['#2c7bb6', '#abd9e9']
    
    # Sort chroms correctly
    df['chrom'] = pd.to_numeric(df['chrom'], errors='coerce')
    df = df.dropna(subset=['chrom'])
    df = df.sort_values(by=['chrom', 'pos'])
    
    # Create continuous x-axis
    df['ind'] = range(len(df))
    df_grouped = df.groupby('chrom')
    
    x_labels = []
    x_labels_pos = []
    
    for num, (name, group) in enumerate(df_grouped):
        ax.scatter(group['ind'], group['-log10(P)'], color=colors_manh[num % len(colors_manh)], s=15, alpha=0.8)
        x_labels.append(name)
        x_labels_pos.append(group['ind'].iloc[0] + (len(group)/2))
    
    # Bonferroni line
    ax.axhline(sig_thresh, color='red', lw=1, linestyle='--', alpha=0.7)
    
    # Set x-ticks
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, rotation=0, fontsize=9)
    ax.set_xlim([0, len(df)])
    
    ax.set_xlabel('Chromosome', fontweight='bold')
    ax.set_ylabel('-log$_{10}$ (p-value)', fontweight='bold')
    ax.set_title('Manhattan Plot of Local Bivariate Genetic Correlations', fontweight='bold', pad=15)
    sns.despine()
    
    plt.savefig(os.path.join(output_dir, 'Manhattan_Plot_LAVA.pdf'))
    plt.savefig(os.path.join(output_dir, 'Manhattan_Plot_LAVA.png'))
    plt.close()
    
    print(f"Finished! Plots saved to {output_dir}")

if __name__ == "__main__":
    import os
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    data_file = os.path.join(BASE_DIR, "LAVA_Local_rg_Robust.csv")
    out_dir = os.path.join(BASE_DIR, "plots")
    create_publication_plots(data_file, out_dir)
