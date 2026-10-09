"""POST-HOC (not pre-registered) expression checks on the K27-LIVER result (DESIGN_V4.md Addendum 6)."""
from pathlib import Path
src = Path("v4k27.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []; LAY")], "k", "exec"))
zz = lambda a: (np.asarray(a, float) - np.nanmean(a)) / np.nanstd(a, ddof=1)
rows = []
for nm, fn in [("GSE17251", "labels_GSE17251.csv"), ("GSE53399", "labels_GSE53399.csv")]:
    L = pd.read_csv(V4 / fn, index_col=0); inL = uniq & pd.Series(sym).isin(L.index).to_numpy()
    expr = pd.Series(np.nan, index=range(len(gh))); expr[inL] = L.expr.reindex(sym[inL]).to_numpy()
    for label, thr in [("all mapped genes, expression as covariate", None), ("top expression quartile", 0.75), ("top expression decile", 0.90)]:
        sel = inL if thr is None else inL & (expr > expr[inL].quantile(thr)).to_numpy()
        idx = np.flatnonzero(sel); st = L.stat.reindex(sym[idx]).to_numpy()
        C = np.column_stack([zz(cov.lpl.to_numpy()[idx]), zz(cov.gc.to_numpy()[idx]), zz(cov.lnt.to_numpy()[idx]), zz(expr.to_numpy()[idx])])
        for frac in (0.05, 0.01):
            y = pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y); r = dict(contrast=f"{nm} top {int(frac*100)}%", subset=label, n=len(idx))
            a0, a1 = SCS["TF-A"][idx], SCS["K27-LIVER"][idx]; df = []
            for _ in range(300):
                a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
                df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
            r.update(auc_tfa=auc_of(a0, y), auc_k27=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5))
            for lab, xx in [("tfa", a0), ("k27", a1)]:
                fm = sm.Logit(y.astype(int), sm.add_constant(np.column_stack([zz(xx), C]))).fit(disp=0); r[f"or_{lab}"] = np.exp(fm.params[1])
            rows.append(r)
o = pd.DataFrame(rows); o.to_csv(V4 / "v4k27_posthoc.tsv", sep="\t", index=False); pd.set_option("display.width", 230); print(o.round(3).to_string(index=False))
