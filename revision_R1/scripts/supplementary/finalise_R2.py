"""Final clean-up of the R2 deliverables: document properties and stray internal text."""
import re
import docx
import openpyxl

FIG = "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/"
AUTHOR = "Vijayachitra Modhukur"
TITLES = {
    "Manuscript_R2_revised.docx": "Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis",
    "Manuscript_R2_changes_highlighted.docx": "Integrative genomic analyses reveal regional genetic sharing and stromal context between uterine fibroids and endometriosis",
    "Response_to_Reviewers_FIG_R2.docx": "Response to the Editor and Reviewers",
    "Supplementary_Figures_FIG_R2.docx": "Supplementary Figures",
}
for f, title in TITLES.items():
    d = docx.Document(FIG + f)
    cp = d.core_properties
    cp.author, cp.last_modified_by, cp.title = AUTHOR, AUTHOR, title
    cp.comments = cp.subject = cp.keywords = cp.category = ""
    d.save(FIG + f)

x = FIG + "Supplementary_Tables_S1_to_S6_F_IG_R2.xlsx"
wb = openpyxl.load_workbook(x)
fix = [
    (r" \(primary; files magma/fibroids_neff\.genes\.out and magma/endo_neff\.genes\.out; source of Figure 3\)", " (primary; source of Figure 3)"),
    (r" \(cell_type_marker_genes\.csv\)", ""),
    (r" Cell-level reannotation code: scripts/reannotation/reannotate_GSE218044\.R in the code repository; per-cell labels: results/reannotation/GSE218044_cell_level_labels\.csv\.gz\.", " The reannotation code and per-cell labels are in the archived code release."),
]
n = 0
for ws in wb:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str):
                v = c.value
                for pat, rep in fix:
                    v = re.sub(pat, rep, v)
                if v != c.value:
                    c.value = v; n += 1
p = wb.properties
p.creator = p.lastModifiedBy = AUTHOR
p.title = "Supplementary Tables S1 to S6"
p.description = p.subject = p.keywords = None
wb.save(x)
print("cell texts cleaned:", n)
