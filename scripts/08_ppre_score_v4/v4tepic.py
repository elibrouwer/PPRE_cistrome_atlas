"""Descriptive head-to-head: original TEPIC-style sum (30-motif run, ChIPseeker same-strand) vs the final score, same evaluation."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
orig = gh.PPRE_score_sum.fillna(0).to_numpy(float); tfa = SC["default"]
plen = cov.lpl.to_numpy(); lnt = cov.lnt.to_numpy()
print("Spearman with promoter length: original", round(stats.spearmanr(orig, plen)[0], 3), "| final", round(stats.spearmanr(tfa, plen)[0], 3))
print("Spearman with #TSS: original", round(stats.spearmanr(orig, lnt)[0], 3), "| final", round(stats.spearmanr(tfa, lnt)[0], 3))
print("Spearman with GC: original", round(stats.spearmanr(orig, cov.gc.to_numpy())[0], 3), "| final", round(stats.spearmanr(tfa, cov.gc.to_numpy())[0], 3))
print("Spearman original vs final:", round(stats.spearmanr(orig, tfa)[0], 3))
rows = []
for nm, idx, y in CON:
    pi, ni = np.flatnonzero(y), np.flatnonzero(~y); a0, a1 = orig[idx], tfa[idx]; df = []
    for _ in range(500):
        a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
        df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
    rows.append(dict(test="hepatocyte " + nm, original=auc_of(a0, y), final=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
pl = np.flatnonzero(ylit); nl = np.flatnonzero(~ylit); dl = []
for _ in range(500):
    a = rng.choice(pl, len(pl)); b = rng.choice(nl, len(nl)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
    dl.append(auc_of(tfa[ii], yy) - auc_of(orig[ii], yy))
rows.append(dict(test="literature genes", original=auc_of(orig, ylit), final=auc_of(tfa, ylit), diff=auc_of(tfa, ylit) - auc_of(orig, ylit), lo=np.percentile(dl, 2.5), hi=np.percentile(dl, 97.5)))
for s_ in SETS:
    ok = (gchr % 2 == 0) if s_ == "PPARG_adipocyte" else None
    o0, o1 = adj(CHIPY[s_], orig, ok), adj(CHIPY[s_], tfa, ok); rows.append(dict(test="ChIP " + s_ + " adjOR", original=o0, final=o1, diff=o1 - o0))
T = pd.DataFrame(rows); T.to_csv(V4 / "v4tepic_headtohead.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(T.round(3).to_string(index=False))
