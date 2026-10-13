"""Sensitivity of the site-count score to the 8 repeat-enriched motifs (DESIGN_V4.md Addendum 16)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
Pn = Pn.assign(q=1.0)
REP = ["MA1556.1", "PPARD_M05796_3.10", "RXRA.H13CORE.0.PS.A", "RXRA_M08962_3.10", "RXRB_M05896_3.10", "RXRB_M10655_3.10", "RXRG_M05828_3.10", "RXRG_M05830_3.10"]
allm = sorted(Hn.motif.unique()); print("motifs in hits:", len(allm), "| REP found:", sum(m in allm for m in REP))
gc = cov.gc.to_numpy(); lnt = cov.lnt.to_numpy(); lpl = cov.lpl.to_numpy()
dec = pd.qcut(pd.Series(lpl), 10, labels=False, duplicates="drop"); gct = pd.qcut(pd.Series(gc), 3, labels=False, duplicates="drop")
def pct_score(motifs):
    n_ = np.bincount(loci_table(1e-5, motifs=motifs).gp, minlength=len(gh)).astype(float)
    return pd.DataFrame({"n": n_, "d": dec, "g": gct}).groupby(["d", "g"]).n.rank(pct=True, method="average").to_numpy() * 100, n_
SCS = {}; NS = {}
for nm, ms in [("all motifs (reference)", None), ("without the 8", [m for m in allm if m not in REP]), ("the 8 alone", REP)]: SCS[nm], NS[nm] = pct_score(ms)
ref = SCS["all motifs (reference)"]; rows = []
for nm, x in SCS.items():
    r = evaluate(x); r["variant"] = nm; r["spearman_vs_ref"] = stats.spearmanr(x, ref)[0]; r["genes_with_sites"] = int((NS[nm] > 0).sum()); r["median_sites"] = float(np.median(NS[nm][NS[nm] > 0])); rows.append(r)
G = pd.DataFrame(rows).set_index("variant"); G.to_csv(V4 / "v4repeat_enriched_results.tsv", sep="\t"); pd.set_option("display.width", 250); print(G.round(3).T.to_string())
hep = [c for c in G.columns if c.startswith("hep")]; chip = [c for c in G.columns if c.startswith("ChIP")]; R0 = "all motifs (reference)"
v = "without the 8"; a = (G.loc[v, hep] - G.loc[R0, hep]).abs().max() <= 0.02; b = abs(G.loc[v, "literature"] - G.loc[R0, "literature"]) <= 0.02; n_ = int(((G.loc[v, chip] / G.loc[R0, chip]) >= 0.90).sum())
print(f"{v}: (a) max hepatocyte |dAUROC| {(G.loc[v, hep] - G.loc[R0, hep]).abs().max():.3f} -> {a} | (b) literature {G.loc[v,'literature'] - G.loc[R0,'literature']:+.3f} -> {b} | (c) ChIP ratio>=0.90 {n_}/6 -> {n_ >= 4} | DROPPABLE: {a and b and n_ >= 4}")
