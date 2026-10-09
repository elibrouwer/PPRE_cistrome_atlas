"""Descriptive: does strength weighting make the percentile score more stable than the plain site count? (not pre-registered)"""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
q_tfa = Pn["q"].copy(); s_ = Pn.s.copy()
lpl = cov.lpl.to_numpy(); gcv = cov.gc.to_numpy(); n = len(gh)
grp = pd.qcut(pd.Series(lpl), 10, labels=False, duplicates="drop").to_numpy() * 3 + pd.qcut(pd.Series(gcv), 3, labels=False, duplicates="drop").to_numpy()
def pct(Rv):
    out = np.empty(n)
    for g_ in np.unique(grp): m = grp == g_; out[m] = pd.Series(Rv[m]).rank(pct=True, method="average").to_numpy() * 100
    return out
W = {"site count": 1.0, "sum of (s-4)": (s_ - 4).clip(0, 5), "strength x decay": q_tfa}
VARS = {"gate 1e-4": dict(gate=1e-4), "gate 1e-6": dict(gate=1e-6), "window -1000..+500": dict(gate=1e-5, dlo=-1000, dhi=500), "window -10000..+3000": dict(gate=1e-5, dlo=-10000, dhi=3000), "25 reliable motifs": dict(gate=1e-5, motifs=RELM)}
rows = []
for wn, w in W.items():
    Pn = Pn.assign(q=w); base = pct(raw_R(loci_table(gate=1e-5), "sum")); eb = evaluate(base); hep = [c for c in eb if c.startswith("hep")]
    for vn, kw in VARS.items():
        x = pct(raw_R(loci_table(**kw), "sum")); e = evaluate(x)
        rows.append(dict(weights=wn, change=vn, spearman=stats.spearmanr(x, base)[0], max_dAUROC=max(abs(e[c] - eb[c]) for c in hep)))
T = pd.DataFrame(rows); T.to_csv(V4 / "v4stab_results.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(T.round(3).pivot(index="change", columns="weights", values="spearman").to_string()); print(T.round(3).pivot(index="change", columns="weights", values="max_dAUROC").to_string())
