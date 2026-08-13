import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

# Set publication style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1,
    "ytick.major.size": 0,
})

def parse_or_ci(val):
    if pd.isna(val) or val == "ns":
        return np.nan, np.nan, np.nan
    try:
        val = str(val).strip()
        or_val, ci_part = val.split(" (")
        ci_part = ci_part.replace(")", "")
        lower, upper = ci_part.split("-")
        return float(or_val), float(lower), float(upper)
    except:
        return np.nan, np.nan, np.nan

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
data_file = os.path.join(repo_root, "results", "Table_2_Comorbidity_Associations_43_Traits.csv")
df = pd.read_csv(data_file)

# Parse columns
df[["Fib_OR", "Fib_LCI", "Fib_UCI"]] = df["Fibroids OR (95% CI)"].apply(parse_or_ci).apply(pd.Series)
df[["Endo_OR", "Endo_LCI", "Endo_UCI"]] = df["Endometriosis OR (95% CI)"].apply(parse_or_ci).apply(pd.Series)

# Define categories and traits EXACTLY as needed, EXCLUDING cancers.
categories = {
    "Reproductive": [
        "PCOS",
        "Adenomyosis",
        "Cervical polyp",
        "Ovarian cyst",
        "Endometrial polyp",
        "Heavy menstrual bleeding",
        "Vulvodynia",
        "Chronic pelvic pain"
    ],
    "Metabolic & Cardiovascular": [
        "NAFLD",
        "Type 2 diabetes",
        "Obesity",
        "Dyslipidaemia",
        "Hypertension"
    ],
    "Pain & Inflammatory": [
        "Interstitial cystitis",
        "SLE",
        "Coeliac disease",
        "IBD",
        "Migraine",
        "Peptic ulcer disease",
        "Rheumatoid arthritis",
        "IBS",
        "Back pain",
        "GERD/reflux",
        "Osteoarthritis"
    ],
    "Neuropsychological": [
        "PTSD",
        "Fibromyalgia",
        "Anxiety",
        "Depression"
    ],
    "Genitourinary & Other": [
        "Thrombocytopenia",
        "Hyperthyroidism",
        "Vitamin D deficiency",
        "Overactive bladder",
        "Urinary incontinence",
        "Recurrent UTI",
        "Osteoporosis",
        "Hypothyroidism"
    ]
}

# Reverse order of categories and traits within categories to plot top-to-bottom
plot_data = []
y_pos = 0
y_ticks = []
y_labels = []

# Gap between categories
category_y_lines = []

for cat in reversed(list(categories.keys())):
    traits = reversed(categories[cat])
    
    first_in_cat = True
    for trait in traits:
        row = df[df["Trait"].str.lower() == trait.lower()]
        if len(row) == 0:
            row = df[df["Trait"].str.contains(trait.split("/")[0], case=False, na=False)]
        
        if len(row) > 0:
            row = row.iloc[0]
            plot_data.append({
                "y": y_pos,
                "trait": trait,
                "category": cat,
                "Fib_OR": row["Fib_OR"], "Fib_LCI": row["Fib_LCI"], "Fib_UCI": row["Fib_UCI"],
                "Endo_OR": row["Endo_OR"], "Endo_LCI": row["Endo_LCI"], "Endo_UCI": row["Endo_UCI"]
            })
            y_ticks.append(y_pos)
            y_labels.append(trait)
            y_pos += 1
            first_in_cat = False
            
    # Draw line above category and add title
    if not first_in_cat:
        category_y_lines.append(y_pos - 0.5)
        plot_data.append({"y": y_pos, "trait": f"__CAT_TITLE__{cat}", "category": cat})
        y_pos += 1

fig, ax = plt.subplots(figsize=(10, 14))

# Plot data
for item in plot_data:
    if item["trait"].startswith("__CAT_TITLE__"):
        cat_title = item["trait"].replace("__CAT_TITLE__", "")
        ax.text(0.5, item["y"], cat_title, ha='left', va='center', fontweight='bold', fontsize=11)
    else:
        y = item["y"]
        # Offset slightly to prevent total overlap
        offset = 0.15
        
        # Endometriosis (Red)
        if not np.isnan(item["Endo_OR"]):
            ax.errorbar(item["Endo_OR"], y + offset, 
                        xerr=[[item["Endo_OR"] - item["Endo_LCI"]], [item["Endo_UCI"] - item["Endo_OR"]]], 
                        fmt='o', color='#e74c3c', capsize=0, elinewidth=1.5, markersize=6)
            
        # Fibroids (Blue)
        if not np.isnan(item["Fib_OR"]):
            ax.errorbar(item["Fib_OR"], y - offset, 
                        xerr=[[item["Fib_OR"] - item["Fib_LCI"]], [item["Fib_UCI"] - item["Fib_OR"]]], 
                        fmt='o', color='#3498db', capsize=0, elinewidth=1.5, markersize=6)

# Aesthetics
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels, fontsize=10)
ax.set_xscale("log")

# Vertical line at 1
ax.axvline(x=1, color='gray', linestyle='--', linewidth=1)

# Horizontal category lines
for y_line in category_y_lines[:-1]:  # Don't draw line at the very top above the last category
    ax.axhline(y=y_line, color='black', linestyle='-', linewidth=1)

ax.set_xlabel("Observational Odds Ratio (OR) [Log Scale]", fontweight='bold', fontsize=12)

# X-axis ticks
ticks = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
ax.set_xticks(ticks)
ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax.set_xlim(0.4, 30)

# Grid
ax.grid(True, axis='x', color='whitesmoke', alpha=0.8, linestyle='-')

# Remove borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

# Custom Legend
import matplotlib.lines as mlines
red_line = mlines.Line2D([], [], color='#e74c3c', marker='o', markersize=7, label='Endometriosis')
blue_line = mlines.Line2D([], [], color='#3498db', marker='o', markersize=7, label='Uterine Fibroids')
ax.legend(handles=[red_line, blue_line], loc='upper center', bbox_to_anchor=(0.5, 1.03), ncol=2, frameon=False, fontsize=11)

plt.tight_layout()

out_dir = os.path.join(repo_root, "figures")
os.makedirs(out_dir, exist_ok=True)
out_file = os.path.join(out_dir, "Figure1_HeadToHead_NoCancers")

plt.savefig(out_file + ".png", dpi=600, bbox_inches='tight')
plt.savefig(out_file + ".pdf", dpi=600, bbox_inches='tight')
print(f"Saved {out_file}.png at 600 DPI.")
