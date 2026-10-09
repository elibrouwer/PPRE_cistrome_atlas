"""Simplest weights under the B2 percentile (DESIGN_V4.md Addendum 14)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
q_tfa = Pn["q"].copy()
gc = cov.gc.to_numpy(); lnt = cov.lnt.to_numpy(); lpl = cov.lpl.to_numpy()
df0 = pd.DataFrame({"lpl": lpl, "gc": gc}); dec = pd.qcut(df0.lpl, 10, labels=False, duplicates="drop"); gct = pd.qcut(df0.gc, 3, labels=False, duplicates="drop")
def pct_score(weight):
    global Pn
    Pn = Pn.assign(q=weight); R_ = raw_R(loci_table(1e-5), "sum"); Pn = Pn.assign(q=q_tfa)
    return pd.DataFrame({"R": R_, "d": dec, "g": gct}).groupby(["d", "g"]).R.rank(pct=True, method="average").to_numpy() * 100
SCS = {"reference B2 (strength x decay)": pct_score(q_tfa), "E2 site count": pct_score(1.0), "E1 sum of -log10 p": pct_score(Pn.s), "E3 sum of (s-4)": pct_score((Pn.s - 4).clip(0, 5))}
ref = SCS["reference B2 (strength x decay)"]; rows = []
for nm, x in SCS.items():
    r = evaluate(x); r["variant"] = nm; r["rho_length"] = stats.spearmanr(x, lpl)[0]; r["rho_nTSS"] = stats.spearmanr(x, lnt)[0]; r["rho_GC"] = stats.spearmanr(x, gc)[0]; r["spearman_vs_ref"] = stats.spearmanr(x, ref)[0]; rows.append(r)
G = pd.DataFrame(rows).set_index("variant"); G.to_csv(V4 / "v4easy2_results.tsv", sep="\t"); pd.set_option("display.width", 260); print(G.round(3).T.to_string())
hep = [c for c in G.columns if c.startswith("hep")]; chip = [c for c in G.columns if c.startswith("ChIP")]; R0 = "reference B2 (strength x decay)"
for v in list(SCS)[1:]:
    a = (G.loc[v, hep] - G.loc[R0, hep]).abs().max() <= 0.02; b = abs(G.loc[v, "literature"] - G.loc[R0, "literature"]) <= 0.02; n_ = int(((G.loc[v, chip] / G.loc[R0, chip]) >= 0.90).sum()); c = n_ >= 4
    mx = max(abs(G.loc[v, "rho_length"]), abs(G.loc[v, "rho_nTSS"]), abs(G.loc[v, "rho_GC"])); d_ = mx <= 0.15
    print(f"{v}: (a) max hepatocyte |dAUROC| {(G.loc[v, hep] - G.loc[R0, hep]).abs().max():.3f} -> {a} | (b) literature {G.loc[v,'literature'] - G.loc[R0,'literature']:+.3f} -> {b} | (c) ChIP ratio>=0.90 {n_}/6 -> {c} | (d) max |rho| {mx:.3f} -> {d_} | EQUIVALENT: {a and b and c and d_}")
