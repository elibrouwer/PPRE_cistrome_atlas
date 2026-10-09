"""H3K27ac layers (DESIGN_V4.md Addendum 6)."""
from pathlib import Path
src = Path("v4_windows.py").read_text(encoding="utf8")
exec(compile(src[:src.index("CFG = {")], "w", "exec"))
cvz = (cov - cov.mean()) / cov.std()
def peaks(path_list):
    d = pd.concat([pd.read_csv(G0 + p, sep="\t", header=None, usecols=[0, 1, 2]) for p in path_list]); d.columns = ["c", "s", "e"]
    d = d[d.c.isin(chrom_ix)].copy(); d["ci"] = d.c.map(chrom_ix); d = d.sort_values(["ci", "s"]).reset_index(drop=True)
    d["pm"] = d.groupby("ci").e.cummax().groupby(d.ci).shift().fillna(-1)
    g = ((d.ci.diff().fillna(1) != 0) | (d.s > d.pm)).cumsum()
    d = d.groupby(g).agg(ci=("ci", "first"), s=("s", "min"), e=("e", "max")).reset_index(drop=True)
    return IntervalSet(d.ci.to_numpy(), d.s.to_numpy(), d.e.to_numpy()), len(d)
LV, nl = peaks(["liver_H3K27ac_ENCFF805YRQ.bed.gz"]); HT, nh_ = peaks(["heart_H3K27ac_ENCFF616PCH.bed.gz"]); UN, nu = peaks(["liver_H3K27ac_ENCFF805YRQ.bed.gz", "heart_H3K27ac_ENCFF616PCH.bed.gz"])
print("merged peaks liver/heart/union:", nl, nh_, nu)
args = (H.cidx.to_numpy(), H.start0.to_numpy(), H.end.to_numpy())
flags = {"K27-LIVER": LV.overlaps(*args), "K27-HEART": HT.overlaps(*args), "K27-UNION": UN.overlaps(*args)}
print({k: round(float(v.mean()), 3) for k, v in flags.items()}, "share of hits in a peak")
def gene_R_from(q, mask):
    p = P[mask]; x = pd.DataFrame({"gene_idx": p.gene_idx.to_numpy(), "locus": p.locus.to_numpy(), "v": np.asarray(q)[mask.to_numpy()].clip(max=0.999999)})
    lg = x.groupby(["gene_idx", "locus"]).v.max().reset_index(); lg["r"] = -np.log1p(-lg.v)
    return lg.groupby("gene_idx").r.sum().reindex(gh.gene_idx).fillna(0).to_numpy()
prom = P.edge <= 0
base = 0.15 * ((P.s - 4) / 5).clip(0, 1).to_numpy() * np.exp(-np.abs(P.d.to_numpy()) / 3000)
SCS = {"TF-A": resid(gene_R_from(base, prom))}
for k, v in flags.items(): SCS[k] = resid(gene_R_from(base * np.where(v[P.hid.to_numpy()], 1.0, 0.5), prom))
print("Spearman with TF-A:", {k: round(stats.spearmanr(v, SCS["TF-A"])[0], 3) for k, v in SCS.items()})
sym = gh.SYMBOL.to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
rows = []; LAY = ["K27-LIVER", "K27-UNION", "K27-HEART"]
for nm, fn in [("GSE17251", "labels_GSE17251.csv"), ("GSE53399", "labels_GSE53399.csv")]:
    L = pd.read_csv(V4 / fn, index_col=0); inL = uniq & pd.Series(sym).isin(L.index).to_numpy()
    expr = pd.Series(np.nan, index=range(len(gh))); expr[inL] = L.expr.reindex(sym[inL]).to_numpy(); exm = (expr > expr[inL].median()).to_numpy()
    idx = np.flatnonzero(exm); st = L.stat.reindex(sym[idx]).to_numpy()
    for frac in (0.05, 0.01):
        y = pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
        for lay in LAY:
            a0, a1 = SCS["TF-A"][idx], SCS[lay][idx]; df = []
            for _ in range(500):
                a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
                df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
            rows.append(dict(test=f"hepatocyte {nm} top {int(frac*100)}%", score=lay, base=auc_of(a0, y), layer=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
sl = pd.read_csv(f"{R}/benchmark/literature_extraction/literature_PPRE_human_prioritised.csv")
sl = sl[sl.tier.isin([1, 2, 3]) & (sl.ref_class == "TSS") & (sl.window_final == "yes")].copy(); sl["symbol"] = sl.Gene_symbol.str.replace(r"\s*\(.*$", "", regex=True).str.strip().str.upper()
ylit = pd.Series(sym).str.upper().isin(set(sl.symbol)).to_numpy()
for lay in LAY: rows.append(dict(test="literature genes AUROC", score=lay, base=auc_of(SCS["TF-A"], ylit), layer=auc_of(SCS[lay], ylit), diff=auc_of(SCS[lay], ylit) - auc_of(SCS["TF-A"], ylit)))
def adj(y, x, ok=None):
    ok = np.ones(len(y), bool) if ok is None else ok
    z = (x[ok] - x[ok].mean()) / x[ok].std(); fm = sm.Logit(y[ok].astype(int), sm.add_constant(np.column_stack([z, cvz.to_numpy()[ok]]))).fit(disp=0); return np.exp(fm.params[1])
for s_ in SETS:
    pk = PK[s_]; lo_ = np.where(th.strand == "+", th.tss - WIN_UP, th.tss - WIN_DOWN); hi_ = np.where(th.strand == "+", th.tss + WIN_DOWN, th.tss + WIN_UP)
    ov = pk.overlaps(th.cidx.to_numpy(), np.maximum(lo_ - 1, 0), hi_); fl = np.zeros(len(gh), bool); fl[np.unique(th.gene_idx.to_numpy()[ov])] = True
    ok = (gchr % 2 == 0) if s_ == "PPARG_adipocyte" else None; o0 = adj(fl, SCS["TF-A"], ok)
    for lay in LAY:
        o1 = adj(fl, SCS[lay], ok); rows.append(dict(test="ChIP " + s_ + " adjOR", score=lay, base=o0, layer=o1, diff=o1 / o0))
o = pd.DataFrame(rows); o.to_csv(V4 / "v4k27_results.tsv", sep="\t", index=False); pd.set_option("display.width", 220); print(o.round(3).to_string(index=False))
for lay in LAY:
    A = o[o.score == lay]; h = A[A.test.str.startswith("hepatocyte")]; a_ = int((h.lo > 0).sum()); lit = A[A.test.str.startswith("literature")].iloc[0]["diff"]; c_ = int((A[A.test.str.startswith("ChIP")]["diff"] >= 1.0).sum())
    print(f"{lay}: (a) lower>0: {a_}/4 (need 3) | (b) literature diff {lit:+.3f} (need >= -0.01) | (c) ChIP OR ratio>=1: {c_}/6 (need 4) -> improves: {a_ >= 3 and lit >= -0.01 and c_ >= 4}")
