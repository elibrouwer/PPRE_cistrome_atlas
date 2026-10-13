"""Empirical ceilings for the hepatocyte AUROC: (1) label reproducibility between the two datasets, (2) ChIP-binding 'oracle' (real PPARG/RXRA peaks near the TSS) as predictor."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
Pn = Pn.assign(q=1.0)
Ls = {nm: pd.read_csv(V4 / fn, index_col=0) for nm, fn in [("GSE17251", "labels_GSE17251.csv"), ("GSE53399", "labels_GSE53399.csv")]}
def labs(nm):
    L = Ls[nm]; inL = uniq & pd.Series(sym).isin(L.index).to_numpy(); expr = pd.Series(np.nan, index=range(len(gh))); expr[inL] = L.expr.reindex(sym[inL]).to_numpy()
    st = pd.Series(np.nan, index=range(len(gh))); st[inL] = L.stat.reindex(sym[inL]).to_numpy(); return inL, (expr > expr[inL].median()).to_numpy(), st.to_numpy()
rows = []
# (1) reproducibility: predict dataset A's top-k% labels with dataset B's statistic, on genes expressed in A and present in B
for a, b in [("GSE17251", "GSE53399"), ("GSE53399", "GSE17251")]:
    inA, exA, stA = labs(a); inB, exB, stB = labs(b); idx = np.flatnonzero(exA & inB & np.isfinite(stB))
    for frac in (0.05, 0.01):
        y = pd.Series(stA[idx]).rank(ascending=False).to_numpy() <= round(frac * len(idx)); rows.append((f"label reproducibility: {b} stat -> {a} top {int(frac*100)}%", auc_of(stB[idx], y), len(idx)))
# (2) ChIP oracle: binding flag (peak within the promoter window) predicting the labels
for nm in Ls:
    inL, ex, st = labs(nm); idx = np.flatnonzero(ex)
    for frac in (0.05, 0.01):
        y = pd.Series(st[idx]).rank(ascending=False).to_numpy() <= round(frac * len(idx))
        for s_ in SETS: rows.append((f"ChIP oracle {s_} -> {nm} top {int(frac*100)}%", auc_of(CHIPY[s_][idx].astype(float), y), len(idx)))
        rows.append((f"site-count score -> {nm} top {int(frac*100)}%", auc_of(SCS_ref[idx], y) if False else np.nan, len(idx)))
R_ = pd.DataFrame(rows, columns=["test", "AUROC", "n_genes"]).dropna(); R_.to_csv(V4 / "v4ceiling_results.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(R_.round(3).to_string())
