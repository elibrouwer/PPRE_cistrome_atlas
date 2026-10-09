"""Simpler score candidates (DESIGN_V4.md Addendum 12)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
def score_with(weight, agg="sum"):
    global Pn
    Pn = Pn.assign(q=weight)
    lg = loci_table(1e-5)
    return resid(raw_R(lg, agg)) if agg != "count" else resid(np.bincount(lg.gp, minlength=len(gh)).astype(float))
tfa = SC["default"]
q_tfa = Pn["q"].copy()
S1 = score_with((Pn.s - 4).clip(0, 5), "sum")
S2 = score_with(1.0, "count")
Pn["q"] = q_tfa
rows = []
for nm, x in [("TF-A", tfa), ("S1 sum of (s-4)", S1), ("S2 site count", S2)]:
    r = evaluate(x); r["variant"] = nm; r["spearman_vs_TFA"] = stats.spearmanr(x, tfa)[0]; rows.append(r)
G = pd.DataFrame(rows).set_index("variant"); G.to_csv(V4 / "v4simple_results.tsv", sep="\t"); pd.set_option("display.width", 250); print(G.round(3).T.to_string())
hep = [c for c in G.columns if c.startswith("hep")]; chip = [c for c in G.columns if c.startswith("ChIP")]
for v in ["S1 sum of (s-4)", "S2 site count"]:
    a = G.loc[v, "spearman_vs_TFA"] >= 0.95; b = (G.loc[v, hep] - G.loc["TF-A", hep]).abs().max() <= 0.01; c = abs(G.loc[v, "literature"] - G.loc["TF-A", "literature"]) <= 0.02
    d_ = int(((G.loc[v, chip] / G.loc["TF-A", chip]) >= 0.95).sum()) >= 4
    print(f"{v}: (a) spearman {G.loc[v,'spearman_vs_TFA']:.3f} -> {a} | (b) max hepatocyte |dAUROC| {(G.loc[v, hep] - G.loc['TF-A', hep]).abs().max():.3f} -> {b} | (c) literature diff {G.loc[v,'literature'] - G.loc['TF-A','literature']:+.3f} -> {c} | (d) ChIP ratio>=0.95: {int(((G.loc[v, chip] / G.loc['TF-A', chip]) >= 0.95).sum())}/6 -> {d_} | EQUIVALENT: {a and b and c and d_}")
