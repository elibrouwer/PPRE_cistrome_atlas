"""Easier size corrections (DESIGN_V4.md Addendum 13)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
lg = loci_table(1e-5); W = raw_R(lg, "sum"); final = SC["default"]
plen = gh.prom_len.clip(lower=1).to_numpy(float); gc = cov.gc.to_numpy(); lnt = cov.lnt.to_numpy(); lpl = cov.lpl.to_numpy()
floor = W[W > 0].min() / 2
A = np.log(np.where(W > 0, W, floor) / (plen / 1000.0))
df = pd.DataFrame({"W": W, "lpl": lpl, "gc": gc})
df["dec"] = pd.qcut(df.lpl, 10, labels=False, duplicates="drop"); df["gct"] = pd.qcut(df.gc, 3, labels=False, duplicates="drop")
B1 = df.groupby("dec").W.rank(pct=True, method="average").to_numpy() * 100
B2 = df.groupby(["dec", "gct"]).W.rank(pct=True, method="average").to_numpy() * 100
SCS = {"final (spline residual)": final, "A density per kb": A, "B1 percentile in length decile": B1, "B2 percentile in length decile x GC tertile": B2}
rows = []
for nm, x in SCS.items():
    r = evaluate(x); r["variant"] = nm; r["rho_length"] = stats.spearmanr(x, lpl)[0]; r["rho_nTSS"] = stats.spearmanr(x, lnt)[0]; r["rho_GC"] = stats.spearmanr(x, gc)[0]; r["spearman_vs_final"] = stats.spearmanr(x, final)[0]; rows.append(r)
G = pd.DataFrame(rows).set_index("variant"); G.to_csv(V4 / "v4easy_results.tsv", sep="\t"); pd.set_option("display.width", 260); print(G.round(3).T.to_string())
hep = [c for c in G.columns if c.startswith("hep")]; chip = [c for c in G.columns if c.startswith("ChIP")]
for v in list(SCS)[1:]:
    a = max(abs(G.loc[v, "rho_length"]), abs(G.loc[v, "rho_nTSS"]), abs(G.loc[v, "rho_GC"])) <= 0.15; b = (G.loc[v, hep] - G.loc["final (spline residual)", hep]).abs().max() <= 0.02
    c = abs(G.loc[v, "literature"] - G.loc["final (spline residual)", "literature"]) <= 0.02; n_ = int(((G.loc[v, chip] / G.loc["final (spline residual)", chip]) >= 0.90).sum()); d_ = n_ >= 4
    print(f"{v}: (a) max |rho| {max(abs(G.loc[v,'rho_length']), abs(G.loc[v,'rho_nTSS']), abs(G.loc[v,'rho_GC'])):.3f} -> {a} | (b) max hepatocyte |dAUROC| {(G.loc[v, hep] - G.loc['final (spline residual)', hep]).abs().max():.3f} -> {b} | (c) literature diff {G.loc[v,'literature'] - G.loc['final (spline residual)','literature']:+.3f} -> {c} | (d) ChIP ratio>=0.90 {n_}/6 -> {d_} | QUALIFIES: {a and b and c and d_}")
