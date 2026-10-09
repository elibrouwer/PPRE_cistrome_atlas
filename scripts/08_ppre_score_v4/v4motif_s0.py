"""Excluding only the RXR half-site models (DESIGN_V4.md Addendum 13)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
SM = pd.read_csv(V4 / "motif_summary_42.csv").set_index("motif"); hit_motifs = set(Hn.motif); ALL = sorted(hit_motifs)
HALF = {m for m in hit_motifs if SM.loc[m, "type"] == "RXR" and SM.loc[m, "length_bp"] <= 14}; S0 = hit_motifs - HALF
print("half-site models removed:", len(HALF), sorted(HALF)); print("S0 motifs kept:", len(S0))
def score_for(m): return resid(raw_R(loci_table(1e-5, motifs=m), "noisy"))
x_all = score_for(None); x0 = score_for(S0); rs = np.random.default_rng(20261013)
# (a) ChIP, all chromosomes
rnd = [score_for(set(rs.choice(ALL, len(S0), replace=False))) for _ in range(30)]; rows = []
for k in SETS:
    o_all, o0 = adj(CHIPY[k], x_all), adj(CHIPY[k], x0); p90 = np.percentile([adj(CHIPY[k], r_) for r_ in rnd], 90)
    rows.append(dict(heldout=k, or_all42=o_all, or_S0=o0, random_p90=p90, beats_all=o0 >= o_all, beats_p90=o0 > p90))
T = pd.DataFrame(rows); T.to_csv(V4 / "v4motif_s0_chip.tsv", sep="\t", index=False); print(T.round(3).to_string(index=False))
ra = int(T.beats_all.sum()) >= 4 and int(T.beats_p90.sum()) >= 3
# (b) hepatocyte
hr = []
for nm, idx, y in CON:
    pi, ni = np.flatnonzero(y), np.flatnonzero(~y); a0, a1 = x_all[idx], x0[idx]; df = []
    for _ in range(500):
        a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
        df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
    hr.append(dict(contrast=nm, auc_all42=auc_of(a0, y), auc_S0=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
H_ = pd.DataFrame(hr); H_.to_csv(V4 / "v4motif_s0_hepatocyte.tsv", sep="\t", index=False); print(H_.round(3).to_string(index=False)); rb = int((H_.lo > 0).sum()) >= 3
# (c) literature
lit0, lit1 = auc_of(x_all, ylit), auc_of(x0, ylit); rc = lit1 >= lit0 - 0.01
# (d) cardiomyocyte GW0742 vs size-matched random subsets
L = pd.read_csv(V4 / "labels_cardio_GW0742.csv").dropna(subset=["stat"]); L = L[L.symbol.isin(sym[uniq])].drop_duplicates("symbol").set_index("symbol")
pos = pd.Series(np.arange(len(gh)), index=sym); ia = pos.reindex(L.index).to_numpy(); ok = ~np.isnan(ia); L = L[ok]; ia = ia[ok].astype(int)
exm = (L.baseMean > L.baseMean.median()).to_numpy(); ci_ = ia[exm]; cst = L.stat.to_numpy()[exm]; cy = pd.Series(cst).rank(ascending=False).to_numpy() <= round(0.05 * len(ci_))
c_all, c0 = auc_of(x_all[ci_], cy), auc_of(x0[ci_], cy)
rc_ = np.array([auc_of(score_for(set(rs.choice(ALL, len(S0), replace=False)))[ci_], cy) for _ in range(100)]); pc = (rc_ < c0).mean() * 100
cls = "worse than size explains" if pc < 5 else ("better than size-matched" if pc > 95 else "explained by size"); rd = cls != "worse than size explains"
print(f"cardiomyocyte GW0742 top 5%: all 42 {c_all:.3f}, S0 {c0:.3f}, random same size mean {rc_.mean():.3f} (5-95%: {np.percentile(rc_,5):.3f}-{np.percentile(rc_,95):.3f}), percentile {pc:.0f} -> {cls}")
print(f"RULE (a) ChIP beats all-42 {int(T.beats_all.sum())}/6, random p90 {int(T.beats_p90.sum())}/6 -> {ra} | (b) hepatocyte lower>0 {int((H_.lo > 0).sum())}/4 -> {rb} | (c) literature {lit1:.3f} vs {lit0:.3f} -> {rc} | (d) cardiomyocyte -> {rd} | ADOPT S0: {ra and rb and rc and rd}")
