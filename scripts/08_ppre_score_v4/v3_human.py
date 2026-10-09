"""v3_human.py -- all-human evaluation of PPRE score v3 (SCORE_V3_DESIGN.md, Addendum 1).
v3-mouse = frozen mouse-trained model (refit identically to v3_run.py); v3-human = same form refit on PPARG adipocyte, odd chromosome indices.
Usage: python v3_human.py   (env PPRE_CACHE_DIR = liver-validation cache; peaks read from HPEAKS)"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
CACHE = Path(os.environ["PPRE_CACHE_DIR"])
HPEAKS = Path(os.environ["HPEAKS"])
sys.path.insert(0, str(HERE.parent / "PPARA_liver_ChIP_validation" / "scripts"))
from common import IntervalSet, WIN_UP, WIN_DOWN   # noqa
OUT = HERE / "tables"
rng = np.random.default_rng(20261008)
t0 = time.time()
def say(m): print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)

def ivset(p): return IntervalSet(p.cidx.to_numpy(), p.start0.to_numpy(), p.end.to_numpy())

def promoter_pairs(sp):
    h = pd.read_parquet(CACHE / sp / "hits.parquet", columns=["cidx", "start0", "end", "score", "pvalue", "hit_gc", "flank_gc", "in_promoter_any"],
                        filters=[("in_promoter_any", "==", True)]).drop(columns="in_promoter_any")
    t = pd.read_parquet(CACHE / sp / "tss.parquet")
    t = t.assign(sgn=np.where(t.strand == "+", 1, -1), key=(t.cidx.to_numpy(np.int64) << 32) + t.tss.to_numpy(np.int64)).sort_values("key").reset_index(drop=True)
    c = ((h.start0 + h.end) // 2).to_numpy(np.int64)
    hk = (h.cidx.to_numpy(np.int64) << 32) + c
    lo = np.searchsorted(t.key.to_numpy(), hk - WIN_UP - 1, "left"); hi = np.searchsorted(t.key.to_numpy(), hk + WIN_UP + 1, "right")
    n = hi - lo
    hid = np.repeat(np.arange(len(h)), n)
    tix = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)])
    tt = t.iloc[tix]
    d = (c[hid] - tt.tss.to_numpy()) * tt.sgn.to_numpy()
    keep = (d >= -WIN_UP) & (d <= WIN_DOWN)
    pr = pd.DataFrame({"hit_id": hid[keep], "gene_idx": tt.gene_idx.to_numpy()[keep], "d": d[keep].astype(float)})
    h["s"] = (-np.log10(h.pvalue.clip(lower=1e-12))).clip(4, 9)
    h["hit_id"] = np.arange(len(h))
    nn = pr.assign(ad=pr.d.abs()).sort_values(["hit_id", "ad"]).drop_duplicates("hit_id").set_index("hit_id")
    hn = h.join(nn[["d"]], on="hit_id").dropna(subset=["d"]).reset_index(drop=True)
    return h, hn, pr

FORM = "y ~ cr(s, df=3) + cr(d, df=5) + cr(flank_gc, df=3) + hit_gc"

class Model:
    def __init__(self, tr):
        self.fit = smf.glm(FORM, data=tr, family=sm.families.Binomial()).fit()
        self.gc = (float(tr.flank_gc.mean()), float(tr.hit_gc.mean()))
    def prob(self, df):
        x = df[["s", "d"]].copy(); x["flank_gc"] = self.gc[0]; x["hit_gc"] = self.gc[1]
        return np.asarray(self.fit.predict(x))

def label(h, P): return P.overlaps(h.cidx.to_numpy(), h.start0.to_numpy(), h.end.to_numpy()).astype(int)

# ---- models
hm, hmn, _ = promoter_pairs("mouse")
hmn["y"] = label(hmn, ivset(pd.read_parquet(CACHE / "peaks" / "MouseLiver_GSE35262.parquet")))
M_mouse = Model(hmn[hmn.cidx % 2 == 1]); say("mouse model fitted")
hh, hhn, prh = promoter_pairs("human"); say(f"human promoter hits {len(hhn):,}")
PK = {p.stem: ivset(pd.read_parquet(p)) for p in HPEAKS.glob("*.parquet")}
PK.pop("PPARG_cardiovascular", None)          # 515 peaks: too few
SETS = ["PPARG_adipocyte", "PPARG_digestive_tract", "PPARG_liver", "RXRA_liver", "RXRA_breast", "RXRA_other_cells"]
hhn["y_train"] = label(hhn, PK["PPARG_adipocyte"])
M_human = Model(hhn.assign(y=hhn.y_train)[hhn.cidx % 2 == 1]); say("human model fitted")
pd.DataFrame({"term": M_human.fit.params.index, "coef": M_human.fit.params.values}).to_csv(OUT / "v3human_model_coefficients.tsv", sep="\t", index=False)
g = pd.DataFrame({"d": np.arange(-3000, 1001, 250.0)}); g["s"] = 6.0; g["flank_gc"] = M_human.gc[0]; g["hit_gc"] = M_human.gc[1]
g["p_human_model"] = M_human.fit.predict(g); g["p_mouse_model"] = M_mouse.fit.predict(g.assign(flank_gc=M_mouse.gc[0], hit_gc=M_mouse.gc[1]))
g[["d", "p_human_model", "p_mouse_model"]].to_csv(OUT / "v3human_distance_kernel.tsv", sep="\t", index=False)

# ---- hit level
def boot(y, scores, blocks, nb=50):
    ub, inv = np.unique(blocks, return_inverse=True)
    idx = [np.where(inv == i)[0] for i in range(len(ub))]
    base = {k: roc_auc_score(y, v) for k, v in scores.items()}
    bs = {k: [] for k in scores}
    for _ in range(nb):
        ii = np.concatenate([idx[j] for j in rng.integers(0, len(ub), len(ub))])
        if y[ii].min() == y[ii].max(): continue
        for k in scores: bs[k].append(roc_auc_score(y[ii], scores[k][ii]))
    return base, {k: (np.percentile(v, 2.5), np.percentile(v, 97.5)) for k, v in bs.items()}

rows = []
def hit_eval(name, P, mask=None):
    d = hhn if mask is None else hhn[mask]
    d = d.reset_index(drop=True); d["y"] = label(d, P)
    pos = d[d.y == 1]; neg = d[d.y == 0]
    if len(pos) < 50: return
    d = pd.concat([pos, neg.sample(min(len(neg), 1_000_000), random_state=1)]).reset_index(drop=True)
    gcf = smf.glm("y ~ cr(flank_gc, df=3) + hit_gc", data=d, family=sm.families.Binomial()).fit()
    sc = {"v3 human-trained": M_human.prob(d), "v3 mouse-trained": M_mouse.prob(d),
          "manuscript FIMO x exp(-d/3000)": (d.score * np.exp(-d.d.abs() / 3000)).to_numpy(), "-log10 p only": d.s.to_numpy(),
          "GC only": np.asarray(gcf.fittedvalues)}
    base, ci = boot(d.y.to_numpy(), sc, (d.cidx.astype(np.int64) * 1000 + d.start0 // 1_000_000).to_numpy())
    for k in sc:
        rows.append(dict(set=name, score=k, auc=base[k], lo=ci[k][0], hi=ci[k][1], n_pos=len(pos)))
    say(f"  {name:34s} " + "  ".join(f"{k.split()[0]}{k.split()[1] if k.startswith('v3') else ''}={base[k]:.3f}" for k in sc))
hit_eval("PPARG_adipocyte (even chr, held out)", PK["PPARG_adipocyte"], hhn.cidx % 2 == 0)
for s in SETS[1:]: hit_eval(s, PK[s])
pd.DataFrame(rows).to_csv(OUT / "v3human_hit_level_auc.tsv", sep="\t", index=False)

# ---- gene level
gh = pd.read_parquet(CACHE / "human" / "genes.parquet"); th = pd.read_parquet(CACHE / "human" / "tss.parquet")
ntss = th.groupby("gene_idx").size().reindex(gh.gene_idx).fillna(1).to_numpy()
gchr = th.groupby("gene_idx").cidx.first().reindex(gh.gene_idx).to_numpy()
cov = pd.DataFrame({"lpl": np.log10(gh.prom_len.clip(lower=1)).to_numpy(), "gc": gh.gc.fillna(gh.gc.median()).to_numpy(), "lnt": np.log(ntss)})

def gene_R(model, gate=1e-5):
    g = hh[hh.pvalue <= gate][["hit_id", "cidx", "start0", "end", "s"]].copy()
    g["c"] = (g.start0 + g.end) // 2
    g = g.sort_values(["cidx", "c"]); g["locus"] = (((g.cidx.diff().fillna(1) != 0) | (g.c.diff().fillna(99) > 10))).cumsum().to_numpy()
    p = prh[prh.hit_id.isin(g.hit_id)].merge(g[["hit_id", "s", "locus"]], on="hit_id")
    p["q"] = model.prob(p)
    lg = p.groupby(["gene_idx", "locus"]).q.max().reset_index()
    lg["r"] = -np.log1p(-lg.q.clip(upper=0.999999))
    return lg.groupby("gene_idx").r.sum().reindex(gh.gene_idx).fillna(0).to_numpy()

def resid(R):
    x = np.log(np.where(R > 0, R, R[R > 0].min() / 2))
    m = smf.ols("x ~ cr(lpl, df=5) + cr(gc, df=5) + cr(lnt, df=4)", data=cov.assign(x=x)).fit()
    z = m.resid.to_numpy(); return (z - z.mean()) / z.std()

SC = {"v3 human-trained": resid(gene_R(M_human)), "v3 mouse-trained": resid(gene_R(M_mouse)),
      "manuscript PPRE_score_sum": gh.PPRE_score_sum.fillna(0).to_numpy(float), "hit count": gh.hit_count.fillna(0).to_numpy(float)}
say("gene scores built")
grows = []
for s in SETS:
    P = PK[s]
    lo = np.where(th.strand == "+", th.tss - WIN_UP, th.tss - WIN_DOWN); hi = np.where(th.strand == "+", th.tss + WIN_DOWN, th.tss + WIN_UP)
    ov = P.overlaps(th.cidx.to_numpy(), np.maximum(lo - 1, 0), hi)
    flag = np.zeros(len(gh), bool); flag[np.unique(th.gene_idx.to_numpy()[ov])] = True
    ok = (gchr % 2 == 0) if s == "PPARG_adipocyte" else np.ones(len(gh), bool)
    for k, v in SC.items():
        vv, yy = v[ok], flag[ok].astype(int)
        z = (vv - vv.mean()) / vv.std()
        Xa = np.column_stack([z] + [(cov[c].to_numpy()[ok] - cov[c].to_numpy()[ok].mean()) / cov[c].to_numpy()[ok].std() for c in cov])
        fa = sm.Logit(yy, sm.add_constant(Xa)).fit(disp=0)
        b, se = fa.params[1], fa.bse[1]
        grows.append(dict(set=s + (" (even chr)" if s == "PPARG_adipocyte" else ""), score=k, n_genes=int(ok.sum()), n_peak=int(yy.sum()), auc=roc_auc_score(yy, vv),
                          adj_or_per_sd=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se)))
    say(f"  {s}: " + "  ".join(f"{r['score'].split()[0]}{r['score'].split()[1]}={r['adj_or_per_sd']:.2f}" for r in grows[-4:]))
pd.DataFrame(grows).to_csv(OUT / "v3human_gene_level.tsv", sep="\t", index=False)
pd.DataFrame({"ENSEMBL": gh.ENSEMBL, "SYMBOL": gh.SYMBOL, "v3_human": SC["v3 human-trained"], "v3_mouse": SC["v3 mouse-trained"]}).to_csv(OUT / "v3human_gene_scores.csv", index=False)
say("done")
