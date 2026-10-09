"""Cardiomyocyte test of the heart H3K27ac layer (DESIGN_V4.md Addendum 7)."""
from pathlib import Path
src = Path("v4k27.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []; LAY")], "k", "exec"))     # SCS: TF-A, K27-LIVER, K27-HEART, K27-UNION
rows = []
for ag in ["WY14643", "Rosiglitazone", "GW0742"]:
    L = pd.read_csv(V4 / f"labels_cardio_{ag}.csv").dropna(subset=["stat"]); L = L[L.symbol.isin(sym[uniq])].drop_duplicates("symbol").set_index("symbol")
    pos = pd.Series(np.arange(len(gh)), index=sym)
    idx_all = pos.reindex(L.index).to_numpy(); ok = ~np.isnan(idx_all); L = L[ok]; idx_all = idx_all[ok].astype(int)
    exm = (L.baseMean > L.baseMean.median()).to_numpy(); idx = idx_all[exm]; st = L.stat.to_numpy()[exm]
    for frac in (0.05, 0.01):
        y = pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
        for lay in ["K27-HEART", "K27-LIVER"]:
            a0, a1 = SCS["TF-A"][idx], SCS[lay][idx]; df = []
            for _ in range(500):
                a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
                df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
            rows.append(dict(agonist=ag, top=f"{int(frac*100)}%", n_expressed=len(idx), score=lay, auc_tfa=auc_of(a0, y), auc_layer=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
o = pd.DataFrame(rows); o.to_csv(V4 / "v4cardio_results.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(o.round(3).to_string(index=False))
P5 = o[(o.top == "5%") & (o.score == "K27-HEART") & o.agonist.isin(["WY14643", "Rosiglitazone"])]
print("RULE: both new-contrast differences positive:", bool((P5["diff"] > 0).all()), "| at least one lower bound > 0:", bool((P5.lo > 0).any()), "->", bool((P5["diff"] > 0).all() and (P5.lo > 0).any()))
