"""v4e test (DESIGN_V4.md Addendum 1): residual fitted within expressed genes vs v3 global residual, evaluated on expressed genes."""
from pathlib import Path
src = Path("v4_windows.py").read_text(encoding="utf8")
exec(compile(src[:src.index("CFG = {")], "w", "exec"))     # hits, pairs, M_human, cov, gh, resid_noTSS, score()
# raw promoter noisy-OR (same code path as score(0, None, False), unresidualised)
p = P[P.edge <= 0]
x = pd.DataFrame({"gene_idx": p.gene_idx.to_numpy(), "locus": p.locus.to_numpy(), "v": p.q.to_numpy().clip(max=0.999999)})
lg = x.groupby(["gene_idx", "locus"]).v.max().reset_index(); lg["r"] = -np.log1p(-lg.v)
Rraw = lg.groupby("gene_idx").r.sum().reindex(gh.gene_idx).fillna(0).to_numpy()
lx = np.log(np.where(Rraw > 0, Rraw, Rraw[Rraw > 0].min() / 2))
def fit_resid(mask):
    d = cov.assign(x=lx)[mask]
    m = smf.ols("x ~ cr(lpl, df=5) + cr(gc, df=5) + cr(lnt, df=4)", data=d).fit()
    pred = m.predict(cov.assign(x=lx)); return (lx - np.asarray(pred))
v3 = resid(Rraw)                                           # global fit (as published v3)
sym = gh.SYMBOL.to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
zz = lambda a: (np.asarray(a, float) - np.nanmean(a)) / np.nanstd(a, ddof=1)
def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
rows = []
for nm, f in [("GSE17251 Wy14643", "labels_GSE17251.csv"), ("GSE53399 GW7647", "labels_GSE53399.csv")]:
    L = pd.read_csv(V4 / f, index_col=0); inL = uniq & pd.Series(sym).isin(L.index).to_numpy()
    expr = pd.Series(np.nan, index=range(len(gh))); expr[inL] = L.expr.reindex(sym[inL]).to_numpy()
    exm = (expr > expr[inL].median()).to_numpy()            # expressed = above median among mapped genes
    v4e = fit_resid(pd.Series(exm).to_numpy())             # fitted on expressed genes only
    idx = np.flatnonzero(exm); st = L.stat.reindex(sym[idx]).to_numpy()
    C = np.column_stack([zz(cov.lpl.to_numpy()[idx]), zz(cov.gc.to_numpy()[idx]), zz(cov.lnt.to_numpy()[idx]), zz(expr.to_numpy()[idx])])
    for frac in (0.05, 0.01):
        y = pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
        a3, a4 = v3[idx], v4e[idx]; df = []
        for _ in range(500):
            a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
            df.append(auc_of(a4[ii], yy) - auc_of(a3[ii], yy))
        r = dict(contrast=f"{nm} top {int(frac*100)}%", n_expressed=len(idx), n_up=int(y.sum()), auc_v3=auc_of(a3, y), auc_v4e=auc_of(a4, y),
                 diff=auc_of(a4, y) - auc_of(a3, y), diff_lo=np.percentile(df, 2.5), diff_hi=np.percentile(df, 97.5))
        for lab, xx in [("v3", a3), ("v4e", a4)]:
            f_ = sm.Logit(y.astype(int), sm.add_constant(np.column_stack([zz(xx), C]))).fit(disp=0); r[f"or_{lab}"] = np.exp(f_.params[1])
        rows.append(r)
o = pd.DataFrame(rows); o["improves"] = o.diff_lo > 0; o.to_csv(V4 / "v4e_results.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(o.round(3).to_string(index=False)); print("contrasts improved:", int(o.improves.sum()), "of", len(o), "-> rule needs >= 3")
