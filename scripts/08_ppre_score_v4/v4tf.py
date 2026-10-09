"""Training-free score test (DESIGN_V4.md Addendum 2)."""
from pathlib import Path
src = Path("v4_windows.py").read_text(encoding="utf8")
exec(compile(src[:src.index("CFG = {")], "w", "exec"))
cvz = (cov - cov.mean()) / cov.std()
def gene_R_from(q, mask):
    p = P[mask]; x = pd.DataFrame({"gene_idx": p.gene_idx.to_numpy(), "locus": p.locus.to_numpy(), "v": np.asarray(q)[mask.to_numpy()].clip(max=0.999999)})
    lg = x.groupby(["gene_idx", "locus"]).v.max().reset_index(); lg["r"] = -np.log1p(-lg.v)
    return lg.groupby("gene_idx").r.sum().reindex(gh.gene_idx).fillna(0).to_numpy()
prom = P.edge <= 0
f = ((P.s - 4) / 5).clip(0, 1).to_numpy()
SCS = {"v3 (learned)": resid(gene_R_from(P.q.to_numpy(), prom)),
       "TF-A (fixed decay 3000)": resid(gene_R_from(0.15 * f * np.exp(-np.abs(P.d.to_numpy()) / 3000), prom)),
       "TF-B (flat)": resid(gene_R_from(0.15 * f, prom))}
print("Spearman with v3:", {k: round(stats.spearmanr(v, SCS["v3 (learned)"])[0], 3) for k, v in SCS.items()})
sym = gh.SYMBOL.to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
zz = lambda a: (np.asarray(a, float) - np.nanmean(a)) / np.nanstd(a, ddof=1)
def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
rows = []
# (1) hepatocyte response, expressed genes
for nm, fn in [("GSE17251", "labels_GSE17251.csv"), ("GSE53399", "labels_GSE53399.csv")]:
    L = pd.read_csv(V4 / fn, index_col=0); inL = uniq & pd.Series(sym).isin(L.index).to_numpy()
    expr = pd.Series(np.nan, index=range(len(gh))); expr[inL] = L.expr.reindex(sym[inL]).to_numpy(); exm = (expr > expr[inL].median()).to_numpy()
    idx = np.flatnonzero(exm); st = L.stat.reindex(sym[idx]).to_numpy()
    for frac in (0.05, 0.01):
        y = pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
        for tf in ("TF-A (fixed decay 3000)", "TF-B (flat)"):
            a3, a4 = SCS["v3 (learned)"][idx], SCS[tf][idx]; df = []
            for _ in range(500):
                a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
                df.append(auc_of(a4[ii], yy) - auc_of(a3[ii], yy))
            rows.append(dict(test=f"hepatocyte {nm} top {int(frac*100)}%", score=tf, v3=auc_of(a3, y), tf=auc_of(a4, y), diff=auc_of(a4, y) - auc_of(a3, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
# (2) literature
sl = pd.read_csv(f"{R}/benchmark/literature_extraction/literature_PPRE_human_prioritised.csv")
sl = sl[sl.tier.isin([1, 2, 3]) & (sl.ref_class == "TSS") & (sl.window_final == "yes")].copy(); sl["symbol"] = sl.Gene_symbol.str.replace(r"\s*\(.*$", "", regex=True).str.strip().str.upper()
ylit = pd.Series(sym).str.upper().isin(set(sl.symbol)).to_numpy()
for tf in ("TF-A (fixed decay 3000)", "TF-B (flat)"):
    rows.append(dict(test="literature genes AUROC", score=tf, v3=auc_of(SCS["v3 (learned)"], ylit), tf=auc_of(SCS[tf], ylit), diff=auc_of(SCS[tf], ylit) - auc_of(SCS["v3 (learned)"], ylit)))
# (3) gene-level binding in 6 human ChIP sets
def adj(y, x, ok=None):
    ok = np.ones(len(y), bool) if ok is None else ok
    z = (x[ok] - x[ok].mean()) / x[ok].std(); fm = sm.Logit(y[ok].astype(int), sm.add_constant(np.column_stack([z, cvz.to_numpy()[ok]]))).fit(disp=0); return np.exp(fm.params[1])
for s_ in SETS:
    pk = PK[s_]; lo_ = np.where(th.strand == "+", th.tss - WIN_UP, th.tss - WIN_DOWN); hi_ = np.where(th.strand == "+", th.tss + WIN_DOWN, th.tss + WIN_UP)
    ov = pk.overlaps(th.cidx.to_numpy(), np.maximum(lo_ - 1, 0), hi_); fl = np.zeros(len(gh), bool); fl[np.unique(th.gene_idx.to_numpy()[ov])] = True
    ok = (gchr % 2 == 0) if s_ == "PPARG_adipocyte" else None
    o3 = adj(fl, SCS["v3 (learned)"], ok)
    for tf in ("TF-A (fixed decay 3000)", "TF-B (flat)"):
        o4 = adj(fl, SCS[tf], ok); rows.append(dict(test="ChIP " + s_ + " adjOR", score=tf, v3=o3, tf=o4, diff=o4 / o3))
o = pd.DataFrame(rows); o.to_csv(V4 / "v4tf_results.tsv", sep="\t", index=False); pd.set_option("display.width", 220); print(o.round(3).to_string(index=False))
A = o[o.score == "TF-A (fixed decay 3000)"]
h = A[A.test.str.startswith("hepatocyte")]; c_a = int((h.lo > -0.01).sum()); lit = A[A.test.str.startswith("literature")].iloc[0]; ch = A[A.test.str.startswith("ChIP")]
c_b = abs(lit["diff"]) <= 0.03; c_c = int((ch["diff"] >= 0.9).sum())
print(f"RULE (a) hepatocyte contrasts with lower bound > -0.01: {c_a} of 4 (need >=3); (b) literature |diff| <= 0.03: {c_b}; (c) ChIP sets with OR ratio >= 0.9: {c_c} of 6 (need >=4)")
print("TF-A equal to v3:", c_a >= 3 and c_b and c_c >= 4)
