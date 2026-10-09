"""Export the final score (number of PPRE sites per promoter, percentile within promoter-length decile x GC tertile) for all genes."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
Pn = Pn.assign(q=1.0)
lg = loci_table(1e-5)                                   # per gene and site (hits p <= 1e-5, window -3000..+1000, sites merged within 10 bp)
n_sites = np.bincount(lg.gp, minlength=len(gh)).astype(float)
df0 = pd.DataFrame({"lpl": cov.lpl.to_numpy(), "gc": cov.gc.to_numpy()})
dec = pd.qcut(df0.lpl, 10, labels=False, duplicates="drop"); gct = pd.qcut(df0.gc, 3, labels=False, duplicates="drop")
pct = pd.DataFrame({"n": n_sites, "d": dec, "g": gct}).groupby(["d", "g"]).n.rank(pct=True, method="average").to_numpy() * 100
out = pd.DataFrame({"ENSEMBL": gh.ENSEMBL, "SYMBOL": gh.SYMBOL, "n_sites": n_sites.astype(int), "promoter_length": gh.prom_len.to_numpy(), "promoter_gc": gh.gc.to_numpy(),
                    "length_decile": dec + 1, "gc_tertile": gct + 1, "PPRE_score_percentile": pct})
out.to_csv(V4 / "PPRE_gene_scores_sitecount.csv", index=False)
print(out.describe().round(2).to_string()); print("genes:", len(out), "| genes without sites:", int((n_sites == 0).sum()), "| groups:", out.groupby(["length_decile", "gc_tertile"]).ngroups)
