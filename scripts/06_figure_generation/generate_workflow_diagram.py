import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

def create_workflow_diagram(output_path):
    # Much wider figure to allow wide, sleek rectangular boxes
    fig, ax = plt.subplots(figsize=(20, 5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Helper to draw rounded boxes
    def draw_box(x, y, width, height, title, subtitle="", text_y_offset=0.0):
        # Create rounded rectangle
        box = patches.FancyBboxPatch(
            (x - width/2, y - height/2), width, height,
            boxstyle="round,pad=0.015,rounding_size=0.03",
            edgecolor='black', facecolor='white', linewidth=1.5, zorder=2
        )
        ax.add_patch(box)
        
        # Add text
        if subtitle:
            ax.text(x, y + 0.02 + text_y_offset, title, ha='center', va='center', 
                    fontsize=12, fontweight='bold', family='sans-serif')
            ax.text(x, y - 0.03 + text_y_offset, subtitle, ha='center', va='center', 
                    fontsize=10, family='sans-serif', linespacing=1.6)
        else:
            ax.text(x, y + text_y_offset, title, ha='center', va='center', 
                    fontsize=12, fontweight='bold', family='sans-serif')

    # Draw Top Box (Data Integration)
    draw_box(0.5, 0.85, 0.85, 0.16, 
             "Data Integration", 
             "Cohorts: UK Biobank + FinnGen R9  |  Multi-omics: GTEx + scRNA-seq Atlases", 
             text_y_offset=0.01)

    # Draw Middle Boxes
    y_mid = 0.5
    h_mid = 0.22  # Shorter boxes (looks more like a landscape rectangle)
    w_mid = 0.18  # Wider boxes
    spacing = 0.19 # Distance between centers
    x_positions = [0.5 + (i - 2) * spacing for i in range(5)]

    steps = [
        ("1. Epidemiology", "Phenotype models\nTemporal ordering"),
        ("2. Architecture", "Global rg (LDSC)\nLocal rg (LAVA)"),
        ("3. Causal Inference", "PheWAS MR\nBidirectional MR"),
        ("4. Prioritization", "SuSiE Coloc\nMAGMA mapping"),
        ("5. Translation", "Proteome-MR\nSingle-cell validation")
    ]

    for x, (title, subtitle) in zip(x_positions, steps):
        draw_box(x, y_mid, w_mid, h_mid, title, subtitle, text_y_offset=0.01)

    # Draw Bottom Box (Prioritized output)
    draw_box(0.5, 0.15, 0.85, 0.16, 
             "Prioritized Output", 
             "134 local rg shared loci  |  Endo \u2192 IBS causality  |  Stromal vs Smooth Muscle divergence", 
             text_y_offset=0.01)

    # Arrows
    # Top to middle (center)
    ax.annotate('', xy=(0.5, y_mid + h_mid/2 + 0.02), xytext=(0.5, 0.85 - 0.15/2 - 0.02),
                arrowprops=dict(arrowstyle='-|>', lw=1.5, color='black'))
    
    # Middle to bottom (center)
    ax.annotate('', xy=(0.5, 0.15 + 0.12/2 + 0.02), xytext=(0.5, y_mid - h_mid/2 - 0.02),
                arrowprops=dict(arrowstyle='-|>', lw=1.5, color='black'))

    # Between middle boxes
    for i in range(4):
        x_start = x_positions[i] + w_mid/2 + 0.01
        x_end = x_positions[i+1] - w_mid/2 - 0.01
        ax.annotate('', xy=(x_end, y_mid), xytext=(x_start, y_mid),
                    arrowprops=dict(arrowstyle='-|>', lw=1.5, color='black'))

    plt.savefig(output_path + '.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(output_path + '.png', bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
    out_dir = os.path.join(repo_root, "figures")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "Figure0_Workflow_Diagram")
    create_workflow_diagram(out_file)
    print(f"Saved workflow diagram to {out_file}.png/pdf")
