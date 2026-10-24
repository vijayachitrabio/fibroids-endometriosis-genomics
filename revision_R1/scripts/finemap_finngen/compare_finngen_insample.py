"""Credible sets per locus from FinnGen R9 in-sample SuSiE, compared with the 1000 Genomes-LD SuSiE,
and credible sets shared between the two disorders. Windows (GRCh38) are the SuSiE lead windows.
Outputs: results/finemap_finngen/FinnGen_insample_vs_our_SuSiE_by_locus.csv,
         results/finemap_finngen/FinnGen_insample_credible_set_members_5loci.csv"""
import pandas as pd
W = pd.read_csv("SuSiE_Locus_Lead_Windows.csv")
old = pd.read_csv("SuSiE_Shared_Loci_Credible_Set_Summary.csv")
rows, mem = [], []
for tr, f, ourtr in (("Uterine fibroids", "CD2_BENIGN_LEIOMYOMA_UTERI", "Uterine_Fibroids"),
                     ("Endometriosis", "N14_ENDOMETRIOSIS", "Endometriosis")):
    c = pd.read_csv(f"finngen_R9_{f}.SUSIE.cred.summary.tsv", sep="\t")
    s = pd.read_csv(f"finngen_R9_{f}.SUSIE.snp.filter.tsv", sep="\t")
    c["chr"] = c.v.str.split(":").str[0]; c["pos"] = c.v.str.split(":").str[1].astype(int)
    for _, w in W.iterrows():
        h = c[(c.chr == str(w.chr)) & (c.pos.between(w.window_start, w.window_end))]
        sm = s.merge(h[["region", "cs"]], on=["region", "cs"]).assign(locus=w.locus, trait=tr)
        mem.append(sm)
        o = old[(old.locus == w.locus) & (old.trait == ourtr)]
        rows.append(dict(locus=w.locus, trait=tr, our_original_cs_L10=len(o), finngen_cs_in_window=len(h),
                         finngen_good_cs=int(h.good_cs.sum()), finngen_cs_sizes=",".join(map(str, h.cs_size)),
                         finngen_lead_variants=",".join(h.v), finngen_lead_p=",".join(f"{p:.1e}" for p in h.p),
                         our_leads_inside_finngen_cs=int(o.lead_bp.isin(sm.position).sum())))
pd.DataFrame(rows).to_csv("FinnGen_insample_vs_our_SuSiE_by_locus.csv", index=False)
pd.concat(mem).to_csv("FinnGen_insample_credible_set_members_5loci.csv", index=False)
