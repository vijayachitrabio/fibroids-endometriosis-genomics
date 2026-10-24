"""Binary gene-window annotation, equivalent to ldsc make_annot.py
(--gene-set-file --gene-coord-file --windowsize 100000 --bimfile):
a SNP is annotated 1 if it lies within [START - w, END + w] of any gene in the set.
Usage: make_annot.py GENES COORD BIM OUT.annot.gz [w]
"""
import gzip, sys
import numpy as np
import pandas as pd
genes_f, coord_f, bim_f, out_f = sys.argv[1:5]
w = int(sys.argv[5]) if len(sys.argv) > 5 else 100000
genes = set(open(genes_f).read().split())
coord = pd.read_csv(coord_f, sep=r"\s+")
coord = coord[coord.GENE.isin(genes)]
bim = pd.read_csv(bim_f, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
chrom = str(bim.CHR.iloc[0])
c = coord[coord.CHR.astype(str) == chrom]
ann = np.zeros(len(bim), dtype=int)
bp = bim.BP.to_numpy()
for s, e in zip(c.START - w, c.END + w):
    ann[(bp >= s) & (bp <= e)] = 1
with gzip.open(out_f, "wt") as fh:
    fh.write("ANNOT\n")
    fh.write("\n".join(map(str, ann)) + "\n")
print(chrom, len(c), "genes", ann.sum(), "SNPs")
