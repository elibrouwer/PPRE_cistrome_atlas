"""Does the liver-ATAC gate gain survive adjustment for expression? (v4 dev, hepatocyte labels, no new scoring)"""
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
rng = np.random.default_rng(20261008)
V4 = "."
S = pd.read_csv("v4_windows_gene_scores.csv"); B = pd.read_csv("gate1_scores_with_atac.csv")[["SYMBOL", "n_tss", "prom_len", "gc"]]
S = S.merge(B, on="SYMBOL"); S = S[S.SYMBOL.map(S.SYMBOL.value_counts()) == 1].set_index("SYMBOL")
z = lambda a: (np.asarray(a, float) - np.nanmean(a)) / np.nanstd(a, ddof=1)
def auc(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
rows = []
for nm, f in [("GSE17251", "labels_GSE17251.csv"), ("GSE53399", "labels_GSE53399.csv")]:
    L = pd.read_csv(f, index_col=0); idx = S.index.intersection(L.index); d = S.loc[idx]; L = L.loc[idx]
    cov = np.column_stack([z(np.log10(d.prom_len)), z(d.gc), z(np.log(d.n_tss)), z(L.expr)])
    for frac in (.05, .01):
        y = (L.stat.rank(ascending=False) <= round(frac * len(d))).to_numpy()
        hi = (L.expr > L.expr.median()).to_numpy()            # expressed half only
        res = {}
        for sc in ["promoter", "promoter + liver ATAC", "liver ATAC at promoter only"]:
            x = d[sc].to_numpy(float)
            f_ = sm.Logit(y.astype(int), sm.add_constant(np.column_stack([z(x), cov]))).fit(disp=0)
            res[sc] = (np.exp(f_.params[1]), np.exp(f_.params[1] - 1.96 * f_.bse[1]), np.exp(f_.params[1] + 1.96 * f_.bse[1]), auc(x[hi], y[hi]))
        # paired bootstrap of AUC difference (gate - promoter) on expressed half
        xa, xb = d["promoter + liver ATAC"].to_numpy(float)[hi], d["promoter"].to_numpy(float)[hi]; yh = y[hi]; pi, ni = np.flatnonzero(yh), np.flatnonzero(~yh); df = []
        for _ in range(500):
            a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
            df.append(auc(xa[ii], yy) - auc(xb[ii], yy))
        for sc, v in res.items(): rows.append(dict(dataset=nm, top=f"{int(frac*100)}%", score=sc, adj_or=v[0], lo=v[1], hi=v[2], auc_expressed_half=v[3]))
        rows.append(dict(dataset=nm, top=f"{int(frac*100)}%", score="DIFF gate - promoter (expressed half AUC)", adj_or=np.mean(df), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
o = pd.DataFrame(rows); o.to_csv("v4_expr_adjust.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(o.round(3).to_string(index=False))
