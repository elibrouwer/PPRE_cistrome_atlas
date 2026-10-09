"""TF-A robustness grid, per-gene intervals, leftover confounds (DESIGN_V4.md Addendum 4)."""
from pathlib import Path
import pyarrow.parquet as pq
src = Path("v4_windows.py").read_text(encoding="utf8")
exec(compile(src[:src.index("# ---- liver ATAC")], "w", "exec"))     # head + gene section (gh, th, cov, resid, PK, SETS, gchr, IntervalSet)
cvz = (cov - cov.mean()) / cov.std()
OUTD = V4
# ---------------- hits within -10000..+3000 of a TSS, gate 1e-4, per chromosome
tt = th.assign(sgn=np.where(th.strand == "+", 1, -1))
gpos = pd.Series(np.arange(len(gh)), index=gh.gene_idx.to_numpy())
HI, PI = [], []; nh = 0
for k in sorted(tt.cidx.unique()):
    tab = pq.read_table(CACHE / "human42" / "hits.parquet", columns=["cidx", "start0", "end", "pvalue", "motif_id"], filters=[("cidx", "==", int(k)), ("pvalue", "<=", 1e-4)]).to_pandas()
    if not len(tab): continue
    tk = tt[tt.cidx == k].sort_values("tss"); ts = tk.tss.to_numpy(np.int64); tg = tk.gene_idx.to_numpy(); sg = tk.sgn.to_numpy()
    c = ((tab.start0 + tab.end) // 2).to_numpy(np.int64)
    lo = np.searchsorted(ts, c - 10000, "left"); hi = np.searchsorted(ts, c + 10000, "right"); n = hi - lo
    hid = np.repeat(np.arange(len(tab)), n); offs = np.cumsum(n) - n; ti = np.arange(n.sum()) - np.repeat(offs, n) + np.repeat(lo, n)
    d = (c[hid] - ts[ti]) * sg[ti]; keep = (d >= -10000) & (d <= 3000)
    pr = pd.DataFrame({"hit_id": hid[keep] + nh, "gene_idx": tg[ti][keep], "d": d[keep].astype(float)})
    hu = np.unique(pr.hit_id.to_numpy() - nh)
    hh_ = pd.DataFrame({"hit_id": hu + nh, "cidx": k, "c": c[hu], "p": tab.pvalue.to_numpy()[hu], "motif": tab.motif_id.astype(str).to_numpy()[hu]})
    HI.append(hh_); PI.append(pr); nh += len(tab)
Hn = pd.concat(HI, ignore_index=True); Pn = pd.concat(PI, ignore_index=True)
Hn["s"] = (-np.log10(Hn.p.clip(lower=1e-12))).clip(4, 9)
Pn["s"] = Hn.set_index("hit_id").s.reindex(Pn.hit_id).to_numpy()
Pn["q"] = 0.15 * ((Pn.s - 4) / 5).clip(0, 1) * np.exp(-np.abs(Pn.d) / 3000)
say(f"hits near TSS {len(Hn):,}; pairs {len(Pn):,}")
rel = pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_reliability.tsv", sep="\t"); RELM = set(rel[rel.reliable].motif)

def loci_table(gate, motifs=None, dlo=-3000, dhi=1000):
    hs = Hn[Hn.p <= gate]
    if motifs is not None: hs = hs[hs.motif.isin(motifs)]
    hs = hs.sort_values(["cidx", "c"]).copy()
    hs["locus"] = (((hs.cidx.diff().fillna(1) != 0) | (hs.c.diff().fillna(99) > 10))).cumsum().to_numpy()
    p = Pn[(Pn.d >= dlo) & (Pn.d <= dhi) & Pn.hit_id.isin(hs.hit_id)].merge(hs[["hit_id", "locus"]], on="hit_id")
    lg = p.groupby(["gene_idx", "locus"]).q.max().reset_index(); lg["gp"] = gpos.reindex(lg.gene_idx).to_numpy(); return lg
def raw_R(lg, agg="noisy"):
    n = len(gh)
    if agg == "noisy": v = -np.log1p(-lg.q.clip(upper=0.999999)); return np.bincount(lg.gp, weights=v, minlength=n)
    if agg == "sum": return np.bincount(lg.gp, weights=lg.q, minlength=n)
    return pd.Series(lg.q.to_numpy()).groupby(lg.gp.to_numpy()).max().reindex(range(n)).fillna(0).to_numpy()
VARS = {"default": dict(gate=1e-5), "gate 1e-4": dict(gate=1e-4), "gate 1e-6": dict(gate=1e-6), "best hit": dict(gate=1e-5, agg="best"), "sum of q": dict(gate=1e-5, agg="sum"),
        "window -1000..+500": dict(gate=1e-5, dlo=-1000, dhi=500), "window -10000..+3000": dict(gate=1e-5, dlo=-10000, dhi=3000), "25 reliable motifs": dict(gate=1e-5, motifs=RELM)}
SC = {}; LG = {}
for nm, kw in VARS.items():
    kw = dict(kw); agg = kw.pop("agg", "noisy"); lg = loci_table(**kw); LG[nm] = lg; SC[nm] = resid(raw_R(lg, agg)); say(f"scored {nm}")

# ---------------- evaluation machinery
sym = gh.SYMBOL.to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
CON = []
for nm, fn in [("GSE17251", "labels_GSE17251.csv"), ("GSE53399", "labels_GSE53399.csv")]:
    L = pd.read_csv(V4 / fn, index_col=0); inL = uniq & pd.Series(sym).isin(L.index).to_numpy()
    expr = pd.Series(np.nan, index=range(len(gh))); expr[inL] = L.expr.reindex(sym[inL]).to_numpy(); exm = (expr > expr[inL].median()).to_numpy()
    idx = np.flatnonzero(exm); st = L.stat.reindex(sym[idx]).to_numpy()
    for frac in (0.05, 0.01): CON.append((f"{nm} top {int(frac*100)}%", idx, pd.Series(st).rank(ascending=False).to_numpy() <= round(frac * len(idx))))
sl = pd.read_csv(f"{R}/benchmark/literature_extraction/literature_PPRE_human_prioritised.csv")
sl = sl[sl.tier.isin([1, 2, 3]) & (sl.ref_class == "TSS") & (sl.window_final == "yes")].copy(); sl["symbol"] = sl.Gene_symbol.str.replace(r"\s*\(.*$", "", regex=True).str.strip().str.upper()
ylit = pd.Series(sym).str.upper().isin(set(sl.symbol)).to_numpy()
CHIPY = {}
for s_ in SETS:
    pk = PK[s_]; lo_ = np.where(th.strand == "+", th.tss - WIN_UP, th.tss - WIN_DOWN); hi_ = np.where(th.strand == "+", th.tss + WIN_DOWN, th.tss + WIN_UP)
    ov = pk.overlaps(th.cidx.to_numpy(), np.maximum(lo_ - 1, 0), hi_); fl = np.zeros(len(gh), bool); fl[np.unique(th.gene_idx.to_numpy()[ov])] = True; CHIPY[s_] = fl
def adj(y, x, ok=None):
    ok = np.ones(len(y), bool) if ok is None else ok
    z = (x[ok] - x[ok].mean()) / x[ok].std(); fm = sm.Logit(y[ok].astype(int), sm.add_constant(np.column_stack([z, cvz.to_numpy()[ok]]))).fit(disp=0); return np.exp(fm.params[1])
def evaluate(x):
    r = {f"hep {nm}": auc_of(x[idx], y) for nm, idx, y in CON}; r["literature"] = auc_of(x, ylit)
    for s_ in SETS: r[f"ChIP {s_}"] = adj(CHIPY[s_], x, (gchr % 2 == 0) if s_ == "PPARG_adipocyte" else None)
    return r
rows = []
for nm, x in SC.items():
    r = evaluate(x); r["variant"] = nm; r["spearman_vs_default"] = stats.spearmanr(x, SC["default"])[0]; rows.append(r)
G = pd.DataFrame(rows).set_index("variant"); G.to_csv(V4 / "v4robust_grid.tsv", sep="\t"); pd.set_option("display.width", 250); print(G.round(3).T.to_string())
hep = [c for c in G.columns if c.startswith("hep")]
for v in G.index[1:]:
    ok = G.loc[v, "spearman_vs_default"] >= 0.8 and (G.loc[v, hep] - G.loc["default", hep]).abs().max() <= 0.02
    print(f"STABLE? {v}: spearman {G.loc[v,'spearman_vs_default']:.3f}, max |dAUROC| {(G.loc[v, hep] - G.loc['default', hep]).abs().max():.3f} -> {ok}")

# ---------------- B. per-gene intervals (Poisson bootstrap over loci)
lg = LG["default"]; v_ = -np.log1p(-lg.q.clip(upper=0.999999)).to_numpy(); gp = lg.gp.to_numpy(); n = len(gh)
R0 = np.bincount(gp, weights=v_, minlength=n); floor = R0[R0 > 0].min() / 2
lx0 = np.log(np.where(R0 > 0, R0, floor)); m0 = smf.ols("x ~ cr(lpl, df=5) + cr(gc, df=5) + cr(lnt, df=4)", data=cov.assign(x=lx0)).fit()
pred = np.asarray(m0.predict(cov.assign(x=lx0))); res0 = lx0 - pred; mu, sd = res0.mean(), res0.std(); z0 = (res0 - mu) / sd
sorted0 = np.sort(z0); NB = 200; pct = np.empty((NB, n), np.float32)
for b in range(NB):
    w = rng.poisson(1.0, len(v_)); Rb = np.bincount(gp, weights=w * v_, minlength=n)
    zb = (np.log(np.where(Rb > 0, Rb, floor)) - pred - mu) / sd; pct[b] = np.searchsorted(sorted0, zb) / n * 100
lo_p, hi_p = np.percentile(pct, 2.5, axis=0), np.percentile(pct, 97.5, axis=0); p0 = np.searchsorted(sorted0, z0) / n * 100
top0 = z0 >= np.sort(z0)[-2000]; instab = (pct >= 100 * (1 - 2000 / n)).mean(axis=0)
out = pd.DataFrame({"ENSEMBL": gh.ENSEMBL, "SYMBOL": gh.SYMBOL, "score_TFA": z0, "percentile": p0, "pct_lo": lo_p, "pct_hi": hi_p})
out.to_csv(V4 / "TFA_gene_scores_with_intervals.csv", index=False)
w_ = hi_p - lo_p
print(f"INTERVALS: median width {np.median(w_):.1f} pct pts; share wider than 30: {(w_ > 30).mean():.3f}; top-2000 genes staying in top 2000 in >=80% draws: {(instab[top0] >= 0.8).mean():.3f}")
print("genes with zero loci:", int((R0 == 0).sum()), "| interval width among genes with >=3 loci:", round(float(np.median(w_[np.bincount(gp, minlength=n) >= 3])), 1))

# ---------------- C. leftover confounds
cp = pd.read_csv(f"{R}/benchmark/results/X_promoter_cpg.csv").set_index("ENSEMBL").reindex(gh.ENSEMBL)
t1 = th.sort_values("tss").groupby("gene_idx").first(); key = (t1.cidx.to_numpy(np.int64) << 32) + t1.tss.to_numpy(np.int64); srt = np.sort(key)
dens = pd.Series(np.searchsorted(srt, key + 50001) - np.searchsorted(srt, key - 50001) - 1, index=t1.index).reindex(gh.gene_idx).to_numpy()
conf = {"CpG obs/exp": cp.cpg_oe.to_numpy(), "CpG island fraction": cp.cgi_fraction.to_numpy(), "gene density (TSS within 50 kb)": dens.astype(float)}
flag = []
for k, v in conf.items():
    ok = np.isfinite(v); rho = stats.spearmanr(SC["default"][ok], v[ok])[0]; print(f"LEFTOVER {k}: rho = {rho:+.3f}");
    if abs(rho) > 0.1: flag.append(k)
print("flagged:", flag)
if flag:
    lx = np.log(np.where(raw_R(LG["default"]) > 0, raw_R(LG["default"]), floor)); d_ = cov.assign(x=lx)
    for k in flag: d_[k.replace(" ", "_").replace("/", "").replace("(", "").replace(")", "")] = pd.Series(np.nan_to_num(conf[k], nan=np.nanmedian(conf[k])), index=d_.index)
    extra = " + ".join(c for c in d_.columns if c not in ("lpl", "gc", "lnt", "x"))
    mp = smf.ols("x ~ cr(lpl, df=5) + cr(gc, df=5) + cr(lnt, df=4) + " + extra, data=d_).fit(); zp = mp.resid.to_numpy(); zp = (zp - zp.mean()) / zp.std()
    e0, e1 = evaluate(SC["default"]), evaluate(zp)
    cmp_ = pd.DataFrame({"TF-A": e0, "TF-A+": e1}); cmp_["diff"] = cmp_["TF-A+"] - cmp_["TF-A"]; print(cmp_.round(3).to_string())
    hepd = cmp_.loc[[c for c in cmp_.index if c.startswith("hep")], "diff"]; litd = cmp_.loc["literature", "diff"]
    print(f"TF-A+ adopted? max hepatocyte drop {hepd.min():+.3f} (limit -0.01), literature diff {litd:+.3f} (limit -0.01) ->", hepd.min() >= -0.01 and litd >= -0.01)
    print("rho after extra correction:", {k: round(stats.spearmanr(zp[np.isfinite(conf[k])], conf[k][np.isfinite(conf[k])])[0], 3) for k in flag})
