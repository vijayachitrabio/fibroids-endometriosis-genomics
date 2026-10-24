# LDSC compatibility patch

The S-LDSC reruns used the CBIIT Python 3 fork of LDSC (v1.0.1), the same tool as the original analysis.
With bitarray >= 3, `ldscore/ldscore.py` line 395 fails because `bitarray.decode()` returns an iterator.
One-line fix applied to a copy of the tool:

    - X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T
    + X = np.array(list(slice.decode(self._bedcode)), dtype="float64").reshape((b, nru)).T

Validation: with this patch, rebuilt annotations and LD scores matched the original files exactly
(max |difference| = 0), and the original S-LDSC results were reproduced to ~13 significant figures.
Reference files: Zenodo record 7768714 (1000G_Phase3_baseline_v1.2_ldscores.tgz,
1000G_Phase3_baselineLD_v2.2_ldscores.tgz; MD5 4a8faf1e… and b261e0ca…). Allele frequencies were
computed with `plink --freq` from the 1000 Genomes Phase 3 EUR PLINK files (489 individuals).
