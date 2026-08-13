import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import os

def generate_radial_map():
    # --- 1. CONFIGURATION & DATA ---
    BASE_PATH = os.path.abspath(os.path.dirname(__file__))
    MR_PATH = os.path.join(BASE_PATH, "MR_PheWAS_Robust_Results.csv")
    LAVA_PATH = os.path.join(BASE_PATH, "LAVA_Local_rg_Robust.csv")
    OUTPUT_PATH = os.path.join(BASE_PATH, "plots/Global_Pleiotropy_Hub_v1.png")
    
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    # Load MR Results
    mr_df = pd.read_csv(MR_PATH)
    # Focus on Endometriosis -> Outcome for this view (or pick best signals)
    mr_df = mr_df[mr_df['Exposure'] == 'Endometriosis'].copy()
    
    # Categorize outcomes
    categories = {
        'irritable_bowel_syndrome': 'Gastrointestinal',
        'endometrial_cancer': 'Malignancy',
        'depression': 'Psychological',
        'anxiety': 'Psychological',
        'pcos': 'Reproductive',
        'ovarian_cyst': 'Reproductive',
        'ovarian_cancer': 'Malignancy',
        'breast_cancer': 'Malignancy',
        'obesity': 'Metabolic',
        'hypertension': 'Cardiovascular',
        'osteoporosis': 'Musculoskeletal',
        'coeliac_disease': 'Autoimmune',
        'migraine': 'Neurological',
        'fibromyalgia': 'Musculoskeletal',
        'endometriosis': 'Reproductive',
        'uterine_fibroids': 'Reproductive'
    }
    mr_df['Category'] = mr_df['Outcome'].map(categories)
    mr_df['-log10P'] = -np.log10(mr_df['IVW_P'] + 1e-300)
    
    # Load LAVA Results (Genetic Hubs)
    lava_df = pd.read_csv(LAVA_PATH)
    lava_df = lava_df[lava_df['sig_FDR'] == True].sort_values('p_val').head(8) # Top 8 genetic hubs
    
    # --- 2. PLOTTING SETUP ---
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(14, 14), facecolor='#0b0e14')
    ax = fig.add_subplot(111, polar=True)
    ax.set_facecolor('#0b0e14')
    
    # Remove grid and spines
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    # Sorting and grouping for the radial spokes
    mr_df = mr_df.sort_values(['Category', '-log10P'])
    
    # Data for the plot
    N = len(mr_df)
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    radii = mr_df['-log10P'].values
    width = (2 * np.pi) / N * 0.7
    
    # Map categories to colors
    cat_colors = {
        'Reproductive': '#ff007f', # Neon Pink
        'Malignancy': '#ff4500',   # OrangeRed
        'Psychological': '#00ffcc', # Cyan
        'Metabolic': '#ffff00',    # Yellow
        'Cardiovascular': '#ff0000', # Red
        'Musculoskeletal': '#99ccff', # Light Blue
        'Autoimmune': '#cc00ff',   # Purple
        'Neurological': '#00ff00', # Green
        'Gastrointestinal': '#ffcc00' # Amber
    }
    colors = [cat_colors.get(cat, '#ffffff') for cat in mr_df['Category']]

    # --- 3. DRAW RADIAL BARS (MR SIGNALS) ---
    bars = ax.bar(theta, radii, width=width, bottom=10, color=colors, alpha=0.9, zorder=3)
    
    # Add text labels for the outer layer
    for t, r, label, cat in zip(theta, radii, mr_df['Outcome'], mr_df['Category']):
        angle_deg = np.rad2deg(t)
        alignment = 'left' if 0 <= angle_deg < 180 else 'right'
        rotation = angle_deg if 0 <= angle_deg < 180 else angle_deg + 180
        
        # Clean label (underscore to space, title case)
        display_name = label.replace('_', ' ').title()
        
        ax.text(t, r + 15, display_name, 
                rotation=rotation, rotation_mode='anchor',
                ha=alignment, va='center', 
                color=cat_colors.get(cat, 'white'), weight='bold', fontsize=9)

    # --- 4. DRAW INNER GENETIC HUB (LAVA) ---
    # Draw a central glowing circle
    circle = patches.Circle((0,0), radius=9, color='#ffffff', alpha=0.1, zorder=1)
    ax.add_artist(circle)
    
    # Plot top genes in the inner ring
    lava_n = len(lava_df)
    lava_theta = np.linspace(0.0, 2 * np.pi, lava_n, endpoint=False)
    
    for i, (t, gene) in enumerate(zip(lava_theta, lava_df['nearest_gene'])):
        ax.text(t, 6, str(gene), ha='center', va='center', color='white', 
                weight='heavy', fontsize=12, alpha=0.9)
        # Glow effect/line to hub
        ax.plot([t, t], [0, 5], color='white', alpha=0.3, lw=1, linestyle='--')

    # --- 5. LEGEND & DECORATIONS ---
    # Center text
    ax.text(0, 0, "Uterine\nFibroids\n&\nEndometriosis\nHub", 
            ha='center', va='center', color='white', weight='black', fontsize=10)

    # Category Legend
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=k, 
                              markerfacecolor=v, markersize=10, linestyle='None') 
                       for k, v in cat_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 1.1),
              title="Systemic Categories", frameon=False, fontsize=10, title_fontsize=11)

    # --- 6. ADD RADIUS GUIDES (Significance Thresholds) ---
    # Bonferroni line approx (-log10(0.05/46) ~ 3)
    ax.add_artist(patches.Circle((0,0), radius=13, fill=False, color='white', linestyle=':', alpha=0.4, lw=0.8))
    ax.text(np.pi/2, 13.5, "Bonferroni Threshold", color='white', alpha=0.5, fontsize=8, ha='center')

    plt.suptitle("GLOBAL PLEIOTROPY ARCHITECTURE", color='white', fontsize=22, weight='black', y=0.95)
    plt.title("Integrating Genetic Correlations (LAVA Hub) and Causal Associations (Mendelian Randomization)", 
              color='gray', fontsize=12, pad=20)

    plt.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight', facecolor='#0b0e14')
    print(f"Radial Hub Map generated successfully: {OUTPUT_PATH}")

from matplotlib.lines import Line2D # Required for legend elements
if __name__ == "__main__":
    generate_radial_map()
