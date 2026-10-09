"""Recomputation for the site-count percentile score: gene table with intervals, robustness grid, leftover confounds, ranked list (DESIGN_V4.md Addendum 15)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
Pn["q"] = 1.0                                   # every merged site counts 1
lpl = cov.lpl.to_numpy(); gcv = cov.gc.to_numpy(); lnt = cov.lnt.to_numpy()
dec = pd.qcut(pd.Series(lpl), 10, labels=False, duplicates="drop").to_numpy(); gct = pd.qcut(pd.Series(gcv), 3, labels=False, duplicates="drop").to_numpy()
grp = dec * 3 + gct; n = len(gh)
def pct(Rv):
    out = np.empty(n)
    for g_ in np.unique(grp):
        m = grp == g_; out[m] = pd.Series(Rv[m]).rank(pct=True, method="average").to_numpy() * 100
    return out
def count_for(**kw):
    lg_ = loci_table(**kw); return np.bincount(lg_.gp, minlength=n).astype(float), lg_
N0, lg0 = count_for(gate=1e-5); S0 = pct(N0)
print("default: genes with 0 sites", int((N0 == 0).sum()), "| median sites", np.median(N0), "| 30 groups sizes min/max", np.bincount(grp).min(), np.bincount(grp).max())
# ---- robustness grid (one change at a time)
VARS = {"default": dict(gate=1e-5), "gate 1e-4": dict(gate=1e-4), "gate 1e-6": dict(gate=1e-6), "window -1000..+500": dict(gate=1e-5, dlo=-1000, dhi=500),
        "window -10000..+3000": dict(gate=1e-5, dlo=-10000, dhi=3000), "25 reliable motifs": dict(gate=1e-5, motifs=RELM)}
SCN = {}
for nm, kw in VARS.items(): SCN[nm] = pct(count_for(**kw)[0])
rows = []
for nm, x in SCN.items():
    r = evaluate(x); r["variant"] = nm; r["spearman_vs_default"] = stats.spearmanr(x, S0)[0]; rows.append(r)
G = pd.DataFrame(rows).set_index("variant"); G.to_csv(V4 / "v4final_robustness.tsv", sep="\t"); pd.set_option("display.width", 250); print(G.round(3).T.to_string())
hep = [c for c in G.columns if c.startswith("hep")]
for v in G.index[1:]:
    ok = G.loc[v, "spearman_vs_default"] >= 0.8 and (G.loc[v, hep] - G.loc["default", hep]).abs().max() <= 0.02
    print(f"STABLE? {v}: spearman {G.loc[v,'spearman_vs_default']:.3f}, max |dAUROC| {(G.loc[v, hep] - G.loc['default', hep]).abs().max():.3f} -> {ok}")
# ---- per-gene intervals: Poisson bootstrap over sites, percentile within the gene's group against the default distribution
gp = lg0.gp.to_numpy(); NB = 200; sorted_by_grp = {g_: np.sort(N0[grp == g_]) for g_ in np.unique(grp)}; draws = np.empty((NB, n), np.float32)
for b in range(NB):
    Rb = np.bincount(gp, weights=rng.poisson(1.0, len(gp)), minlength=n); pb = np.empty(n)
    for g_, srt in sorted_by_grp.items():
        m = grp == g_; lo_ = np.searchsorted(srt, Rb[m], "left"); hi_ = np.searchsorted(srt, Rb[m], "right"); pb[m] = (lo_ + hi_) / 2 / len(srt) * 100
    draws[b] = pb
plo, phi = np.percentile(draws, 2.5, axis=0), np.percentile(draws, 97.5, axis=0); w_ = phi - plo
order = np.lexsort((-N0, -S0)); top = np.zeros(n, bool); top[order[:2000]] = True; thr = S0[order[1999]]
stay = (draws >= thr).mean(axis=0)
print(f"INTERVALS: median width {np.median(w_):.1f} pct pts; share wider than 30: {(w_ > 30).mean():.3f}; top-2000 genes staying above the cutoff percentile ({thr:.1f}) in >=80% draws: {(stay[top] >= 0.8).mean():.3f}")
print("median width among genes with >=3 sites:", round(float(np.median(w_[N0 >= 3])), 1), "| genes with >=3 sites:", int((N0 >= 3).sum()))
out = pd.DataFrame({"ENSEMBL": gh.ENSEMBL, "SYMBOL": gh.SYMBOL, "n_sites": N0, "percentile": S0, "pct_lo": plo, "pct_hi": phi, "group": grp})
out.to_csv(V4 / "PPRE_sitecount_gene_scores.csv", index=False)
# ---- leftover confounds
cp = pd.read_csv(f"{R}/benchmark/results/X_promoter_cpg.csv").set_index("ENSEMBL").reindex(gh.ENSEMBL)
t1 = th.sort_values("tss").groupby("gene_idx").first(); key = (t1.cidx.to_numpy(np.int64) << 32) + t1.tss.to_numpy(np.int64); srt = np.sort(key)
dens = pd.Series(np.searchsorted(srt, key + 50001) - np.searchsorted(srt, key - 50001) - 1, index=t1.index).reindex(gh.gene_idx).to_numpy().astype(float)
for k, v in {"promoter length": lpl, "ln #TSS": lnt, "promoter GC": gcv, "CpG obs/exp": cp.cpg_oe.to_numpy(), "CpG island fraction": cp.cgi_fraction.to_numpy(), "gene density (50 kb)": dens}.items():
    ok = np.isfinite(v); print(f"LEFTOVER rho {k}: {stats.spearmanr(S0[ok], v[ok])[0]:+.3f}")
# ---- ranked list for the snRNA-seq runs: same gene set as the existing v3 list
base = pd.read_csv(f"{V3}/MyDRE_package/PPRE_v3_ranked.bed", sep="\t")
lst = base[["ensembl_genes", "symbol"]].merge(out[["ENSEMBL", "percentile"]].rename(columns={"ENSEMBL": "ensembl_genes", "percentile": "score"}), on="ensembl_genes", how="left")
print("ranked list genes:", len(lst), "| missing score:", int(lst.score.isna().sum()))
lst.dropna().to_csv(f"{V3}/MyDRE_package/PPRE_sitecount_ranked.bed", sep="\t", index=False)
