"""Cardiomyocyte GW0742 confirmation of restricted motif sets (DESIGN_V4.md Addendum 11, 'confirmation of the PPAR-centric motif set').
NOTE: this file was RECONSTRUCTED on 2026-10-08 after the original was overwritten by accident by another script of the same name. It follows the protocol
written in Addendum 11 (GSE160987 GW0742 vs Control, expressed genes = baseMean above the median, up = top 5% / 1% by Wald statistic, paired AUROC
difference vs all 42 motifs, 500 bootstrap resamples over genes) and writes to v4motif3_cardio_reproduced.tsv so the original table v4motif3_cardio.tsv is untouched."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
SM = pd.read_csv(V4 / "motif_summary_42.csv").set_index("motif")
RELM_ALL = set(SM.index[SM.reliable_25]) & set(Hn.motif)
rel = SM.loc[sorted(RELM_ALL)]
SETS_C = {"S2 PPAR-centric": set(rel.index[rel.type.isin(["PPAR", "heterodimer"])]),
          "S1 no RXR half-sites": set(rel.index[~((rel.type == "RXR") & (rel.length_bp <= 14))]),
          "reliable-25": set(rel.index)}
for k, v in SETS_C.items(): print(k, len(v))
def score_for(motifs): return resid(raw_R(loci_table(1e-5, motifs=motifs), "noisy"))
x_all = score_for(None); xs = {k: score_for(v) for k, v in SETS_C.items()}
sym = gh.SYMBOL.to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
pos = pd.Series(np.arange(len(gh)), index=sym)
L = pd.read_csv(V4 / "labels_cardio_GW0742.csv").dropna(subset=["stat"]); L = L[L.symbol.isin(sym[uniq])].drop_duplicates("symbol").set_index("symbol")
idx_all = pos.reindex(L.index).to_numpy(); ok = ~np.isnan(idx_all); L = L[ok]; idx_all = idx_all[ok].astype(int)
exm = (L.baseMean > L.baseMean.median()).to_numpy(); idx = idx_all[exm]; st = L.stat.to_numpy()[exm]; rows = []
for frac in (0.05, 0.01):
    y = pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
    for k, xv in xs.items():
        a0, a1 = x_all[idx], xv[idx]; df = []
        for _ in range(500):
            a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
            df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
        rows.append(dict(top=f"{int(frac*100)}%", set=k, n_expressed=len(idx), auc_all42=auc_of(a0, y), auc_set=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
o = pd.DataFrame(rows); o.to_csv(V4 / "v4motif3_cardio_reproduced.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(o.round(3).to_string(index=False))
orig = V4 / "v4motif3_cardio.tsv"
if orig.exists(): print("ORIGINAL TABLE:"); print(pd.read_csv(orig, sep="\t").round(3).to_string(index=False))
