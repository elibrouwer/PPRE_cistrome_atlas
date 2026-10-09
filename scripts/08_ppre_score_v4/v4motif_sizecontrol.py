"""Size-matched random-subset control for the cardiomyocyte and hepatocyte motif results (DESIGN_V4.md Addendum 12).
Also contains the cardiomyocyte confirmation of Addendum 11 (the original v4motif3.py was replaced on disk by another analysis)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
SM = pd.read_csv(V4 / "motif_summary_42.csv").set_index("motif"); hit_motifs = set(Hn.motif); ALL = sorted(hit_motifs)
RELA = set(pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_reliability.tsv", sep="\t").query("reliable").motif) & hit_motifs
S2 = {m for m in RELA if SM.loc[m, "type"] in ("PPAR", "heterodimer")}; S1 = {m for m in RELA if not (SM.loc[m, "type"] == "RXR" and SM.loc[m, "length_bp"] <= 14)}
SEL = {"S2 PPAR-centric (9)": S2, "S1 no RXR half-sites (20)": S1, "reliable-25": RELA}
def score_for(m): return resid(raw_R(loci_table(1e-5, motifs=m), "noisy"))
L = pd.read_csv(V4 / "labels_cardio_GW0742.csv").dropna(subset=["stat"]); L = L[L.symbol.isin(sym[uniq])].drop_duplicates("symbol").set_index("symbol")
pos = pd.Series(np.arange(len(gh)), index=sym); ia = pos.reindex(L.index).to_numpy(); ok = ~np.isnan(ia); L = L[ok]; ia = ia[ok].astype(int)
exm = (L.baseMean > L.baseMean.median()).to_numpy(); cidx_ = ia[exm]; cst = L.stat.to_numpy()[exm]
cy = pd.Series(cst).rank(ascending=False).to_numpy() <= round(0.05 * len(cidx_))
def metrics(x):
    c = auc_of(x[cidx_], cy); h = np.mean([auc_of(x[idx], y) for nm, idx, y in CON]); return c, h
x_all = score_for(None); base = metrics(x_all); print("all 42: cardio %.3f, hepatocyte mean %.3f" % base)
rs = np.random.default_rng(20261012); rows = []
for nm, mset in SEL.items():
    n = len(mset); c, h = metrics(score_for(mset))
    rnd = np.array([metrics(score_for(set(rs.choice(ALL, n, replace=False)))) for _ in range(100)])
    pc, ph = (rnd[:, 0] < c).mean() * 100, (rnd[:, 1] < h).mean() * 100
    cls = lambda p: "worse than size explains" if p < 5 else ("better than size-matched" if p > 95 else "explained by size")
    rows.append(dict(set=nm, n=n, cardio_set=c, cardio_random_mean=rnd[:, 0].mean(), cardio_random_p5=np.percentile(rnd[:, 0], 5), cardio_random_p95=np.percentile(rnd[:, 0], 95), cardio_percentile=pc, cardio_class=cls(pc),
                     hep_set=h, hep_random_mean=rnd[:, 1].mean(), hep_random_p5=np.percentile(rnd[:, 1], 5), hep_random_p95=np.percentile(rnd[:, 1], 95), hep_percentile=ph, hep_class=cls(ph)))
    print(f"{nm}: cardio {c:.3f} (random mean {rnd[:,0].mean():.3f}, 5-95% {np.percentile(rnd[:,0],5):.3f}-{np.percentile(rnd[:,0],95):.3f}, pct {pc:.0f}) -> {cls(pc)} | hepatocyte mean {h:.3f} (random mean {rnd[:,1].mean():.3f}, pct {ph:.0f}) -> {cls(ph)}")
o = pd.DataFrame(rows); o.to_csv(V4 / "v4motif_sizecontrol.tsv", sep="\t", index=False)
