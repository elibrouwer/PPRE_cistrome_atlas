"""v4 development: wider gene-linking window and liver-ATAC gating, tested on human hepatocyte agonist responses (GSE17251, GSE53399).
Variants fixed before scoring:
 window: promoter (v3 site model, d in -3000..+1000) | 10kb (decay length 3 kb beyond promoter edge, reach 7 kb) | 50kb (decay 10 kb, reach 47 kb)
 gate:   none | liver ATAC (ENCFF488BRH): hit multiplier 1 inside a peak, 0.25 outside
 Outside the promoter, site probability = v3 model at fixed d=-500 (strength only) times exp(-edge distance / decay) times gate.
Aggregation = v3 (loci within 10 bp merged, noisy-OR, length+GC residual without TSS term)."""
import os, sys, gzip
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
from scipy import stats
R = "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code"
V3 = Path(R) / "results/PPRE_score_v3"; V4 = Path(R) / "results/PPRE_score_v4"
G0 = "C:/Users/brouw/AppData/Local/Temp/claude/C--Users-brouw-OneDrive---Universiteit-Utrecht-Master-Bioinformatics-Minor-internship-R-code/f7175609-6765-498e-8602-f582703f131b/scratchpad/gate0/"
__file__ = str(V3 / "v3_human.py")
src = (V3 / "v3_human.py").read_text(encoding="utf8")
exec(compile(src[:src.index("# ---- hit level")], "head", "exec"))
exec(compile(src[src.index("# ---- gene level"):src.index("SC = {")], "gene", "exec"))
OUT = V4
rng = np.random.default_rng(20261008)

# ---- liver ATAC
chrom_ix = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
la = pd.read_csv(G0 + "liver_ATAC_ENCFF488BRH.bed.gz", sep="\t", header=None, usecols=[0, 1, 2]); la.columns = ["c", "s", "e"]
la = la[la.c.isin(chrom_ix)]
la["ci"] = la.c.map(chrom_ix); la = la.sort_values(["ci", "s"]).reset_index(drop=True)
la["pm"] = la.groupby("ci").e.cummax().groupby(la.ci).shift().fillna(-1)
grp = ((la.ci.diff().fillna(1) != 0) | (la.s >= la.pm)).cumsum()   # touching/overlapping intervals merge (>= keeps abutting separate only if disjoint)
grp = ((la.ci.diff().fillna(1) != 0) | (la.s > la.pm)).cumsum()
la = la.groupby(grp).agg(ci=("ci", "first"), s=("s", "min"), e=("e", "max")).reset_index(drop=True)
assert ((la.ci.diff() != 0) | (la.s > la.e.shift())).all()
LIV = IntervalSet(la.ci.to_numpy(), la.s.to_numpy(), la.e.to_numpy())
say(f"liver ATAC peaks {len(la):,}")

# ---- hits (42 motifs, gate 1e-5) paired with all TSS within 50 kb
H = pd.read_parquet(CACHE / "human42" / "hits.parquet", columns=["cidx", "start0", "end", "pvalue"])
H = H[H.pvalue <= 1e-5].reset_index(drop=True)
H["s"] = (-np.log10(H.pvalue.clip(lower=1e-12))).clip(4, 9); H["hit_id"] = np.arange(len(H))
H["c"] = (H.start0 + H.end) // 2
H["open"] = LIV.overlaps(H.cidx.to_numpy(), H.start0.to_numpy(), H.end.to_numpy())
H = H.sort_values(["cidx", "c"]).reset_index(drop=True)
H["locus"] = (((H.cidx.diff().fillna(1) != 0) | (H.c.diff().fillna(99) > 10))).cumsum().to_numpy()
say(f"hits {len(H):,}, in liver ATAC {H.open.mean():.3f}")
t = th.assign(sgn=np.where(th.strand == "+", 1, -1), key=(th.cidx.to_numpy(np.int64) << 32) + th.tss.to_numpy(np.int64)).sort_values("key").reset_index(drop=True)
W = 50000
c = H.c.to_numpy(np.int64); hk = (H.cidx.to_numpy(np.int64) << 32) + c
lo = np.searchsorted(t.key.to_numpy(), hk - W - 1, "left"); hi = np.searchsorted(t.key.to_numpy(), hk + W + 1, "right"); n = hi - lo
offs = np.cumsum(n) - n
hid = np.repeat(np.arange(len(H)), n); tix = np.arange(n.sum()) - np.repeat(offs, n) + np.repeat(lo, n)
d = (c[hid] - t.tss.to_numpy()[tix]) * t.sgn.to_numpy()[tix]
k = np.abs(d) <= W
P = pd.DataFrame({"hid": hid[k], "gene_idx": t.gene_idx.to_numpy()[tix][k], "d": d[k].astype(float)})
del hid, tix, d, k
P["s"] = H.s.to_numpy()[P.hid]; P["locus"] = H.locus.to_numpy()[P.hid]; P["open"] = H.open.to_numpy()[P.hid]
inprom = (P.d >= -3000) & (P.d <= 1000)
P["edge"] = np.where(inprom, 0.0, np.maximum(-3000 - P.d, P.d - 1000))
q = np.empty(len(P)); q[inprom.to_numpy()] = M_human.prob(P[inprom][["s", "d"]])
q[~inprom.to_numpy()] = M_human.prob(P[~inprom][["s"]].assign(d=-500.0))
P["q"] = q
say(f"pairs {len(P):,}, promoter pairs {int(inprom.sum()):,}")

def resid_noTSS(Rv):
    x = np.log(np.where(Rv > 0, Rv, Rv[Rv > 0].min() / 2))
    m = smf.ols("x ~ cr(lpl, df=5) + cr(gc, df=5)", data=cov.assign(x=x)).fit(); z = m.resid.to_numpy(); return (z - z.mean()) / z.std()

def score(reach, decay, gate):
    p = P[P.edge <= reach]
    w = np.exp(-p.edge.to_numpy() / decay) if decay else np.ones(len(p))
    m = np.where(p.open.to_numpy(), 1.0, 0.25) if gate else 1.0
    x = pd.DataFrame({"gene_idx": p.gene_idx.to_numpy(), "locus": p.locus.to_numpy(), "v": (p.q.to_numpy() * w * m).clip(max=0.999999)})
    lg = x.groupby(["gene_idx", "locus"]).v.max().reset_index(); lg["r"] = -np.log1p(-lg.v)
    return resid_noTSS(lg.groupby("gene_idx").r.sum().reindex(gh.gene_idx).fillna(0).to_numpy())

CFG = {"promoter": (0, None), "10 kb": (7000, 3000), "50 kb": (47000, 10000)}
X = {}
for wn, (reach, dec) in CFG.items():
    for gate in (False, True):
        X[f"{wn}{' + liver ATAC' if gate else ''}"] = score(reach, dec, gate); say(f"scored {wn} gate={gate}")
# gene-level baselines
prom_lo = np.where(th.strand == "+", th.tss - 3000, th.tss - 1000); prom_hi = np.where(th.strand == "+", th.tss + 1000, th.tss + 3000)
ov = LIV.overlaps(th.cidx.to_numpy(), np.maximum(prom_lo - 1, 0), prom_hi); fl = np.zeros(len(gh)); fl[np.unique(th.gene_idx.to_numpy()[ov])] = 1
X["liver ATAC at promoter only"] = fl
X["GC only"] = gh.gc.fillna(gh.gc.median()).to_numpy(); X["promoter length only"] = np.log10(gh.prom_len.clip(lower=1)).to_numpy()
X["original sum"] = gh.PPRE_score_sum.fillna(0).to_numpy(float)

# ---- evaluation on hepatocyte response labels + literature genes
z = lambda a: (np.asarray(a, float) - np.nanmean(a)) / np.nanstd(a, ddof=1)
sym = gh.SYMBOL.to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
def aucci(x, y, nb=500):
    pi, ni = np.flatnonzero(y), np.flatnonzero(~y); bs = []
    for _ in range(nb):
        a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); bs.append(auc_of(np.r_[x[a], x[b]], np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]))
    return auc_of(x, y), np.percentile(bs, 2.5), np.percentile(bs, 97.5)
sets = {}
for nm, f in [("GSE17251 Wy14643", "labels_GSE17251.csv"), ("GSE53399 GW7647", "labels_GSE53399.csv")]:
    L = pd.read_csv(V4 / f, index_col=0); idx = np.flatnonzero(uniq & pd.Series(sym).isin(L.index).to_numpy())
    st = L.stat.reindex(sym[idx]).to_numpy(); sets[nm] = (idx, st)
sl = pd.read_csv(f"{R}/benchmark/literature_extraction/literature_PPRE_human_prioritised.csv")
sl = sl[sl.tier.isin([1, 2, 3]) & (sl.ref_class == "TSS") & (sl.window_final == "yes")].copy(); sl["symbol"] = sl.Gene_symbol.str.replace(r"\s*\(.*$", "", regex=True).str.strip().str.upper()
ylit = pd.Series(sym).str.upper().isin(set(sl.symbol)).to_numpy()
rows = []
for name, x in X.items():
    for sn, (idx, st) in sets.items():
        for frac in (0.05, 0.01):
            kk = int(round(frac * len(idx))); y = pd.Series(st).rank(ascending=False).to_numpy() <= kk
            a, l, h_ = aucci(np.asarray(x, float)[idx], y); rows.append(dict(test=f"{sn} top {int(frac*100)}%", score=name, auc=a, lo=l, hi=h_))
    a, l, h_ = aucci(np.asarray(x, float), ylit); rows.append(dict(test="literature genes tiers 1-3", score=name, auc=a, lo=l, hi=h_))
o = pd.DataFrame(rows); o.to_csv(OUT / "v4_windows_results.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(o.assign(r=lambda d: d.auc.round(3).astype(str) + " (" + d.lo.round(3).astype(str) + "-" + d.hi.round(3).astype(str) + ")").pivot(index="score", columns="test", values="r").to_string())
pd.DataFrame({"SYMBOL": gh.SYMBOL, **{k: v for k, v in X.items()}}).to_csv(OUT / "v4_windows_gene_scores.csv", index=False)
