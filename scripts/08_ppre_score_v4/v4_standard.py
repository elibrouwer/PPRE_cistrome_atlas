"""Motif prioritisation re-tested with field-standard benchmarking (DESIGN_V4.md Addendum 15)."""
from pathlib import Path
from sklearn.metrics import average_precision_score
from statsmodels.stats.contingency_tables import StratifiedTable
src = Path("v4robust.py").read_text(encoding="utf8")
src = src.replace('columns=["cidx", "start0", "end", "pvalue", "motif_id"]', 'columns=["cidx", "start0", "end", "pvalue", "motif_id", "flank_gc"]')
src = src.replace('"motif": tab.motif_id.astype(str).to_numpy()[hu]}', '"motif": tab.motif_id.astype(str).to_numpy()[hu], "fgc": tab.flank_gc.to_numpy()[hu]}')
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
INV = pd.read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_human_only.csv", encoding="utf-8-sig").set_index("motif_id")
SM = pd.read_csv(V4 / "motif_summary_42.csv").set_index("motif")
MT = pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_gc_band_auc_by_set.tsv", sep="\t")
A = MT.groupby(["motif", "set"]).auc.mean().unstack(); O = MT[MT.split == "odd"].set_index(["motif", "set"])
hit_motifs = set(Hn.motif); ALL = sorted(hit_motifs)
def select_excluding(k):
    cols = [c for c in A.columns if c != k]; sub = A[cols]; nh = O.n_hits.unstack()[cols]; npz = O.n_pos.unstack()[cols]
    ok = (sub.median(axis=1) >= 0.55) & ((sub > 0.5).sum(axis=1) >= 4) & (nh.median(axis=1) >= 2000) & (npz.median(axis=1) >= 50)
    return set(ok[ok].index) & hit_motifs
def s2(m): return {x for x in m if SM.loc[x, "type"] in ("PPAR", "heterodimer")}
def s1(m): return {x for x in m if not (SM.loc[x, "type"] == "RXR" and SM.loc[x, "length_bp"] <= 14)}
RELM_ALL = set(SM.index[SM.reliable_25]) & hit_motifs
dr1 = set(INV.index[INV.architecture == "DR1"]) & hit_motifs; dr1q = {m for m in dr1 if str(INV.loc[m, "sparse_flag"]) != "SPARSE"}
sub_ = INV.loc[sorted(dr1q)].assign(wid=INV.loc[sorted(dr1q), "width"]).sort_values(["ic_bits", "wid"], ascending=False); nr = set(sub_.groupby("cluster_r098").head(1).index)
s0 = {m for m in hit_motifs if not (SM.loc[m, "type"] == "RXR" and SM.loc[m, "length_bp"] <= 14)}
# sets: name -> (function giving the set for held-out k or None for fixed set)
FIXED = {"all 42": set(ALL), "S0 remove RXR half-sites only (28)": s0, "DR1 (22)": dr1, "DR1-Q (20)": dr1q, "DR1-NR (5)": nr}
LODO = {"reliable-25": lambda k: select_excluding(k), "PPAR-centric (9)": lambda k: s2(select_excluding(k)), "no RXR half-sites (20)": lambda k: s1(select_excluding(k))}

# ---------------- (1) hit level: Mantel-Haenszel OR retained vs removed, flank-GC strata
Hn["band"] = pd.qcut(Hn.fgc, 10, labels=False, duplicates="drop")
hs = Hn[Hn.p <= 1e-5].copy()
LAB = {s: PK[s].overlaps(hs.cidx.to_numpy(), hs.c.to_numpy(), hs.c.to_numpy() + 1) for s in SETS}
def mh(retained_mask, y):
    tabs = []
    for b, g in hs.groupby("band"):
        r = retained_mask[g.index]; yy = y[g.index]
        a, bb, c, d = (yy & r).sum(), ((~yy) & r).sum(), (yy & ~r).sum(), ((~yy) & ~r).sum()
        if min(a + bb, c + d) > 0 and (a + c) > 0: tabs.append(np.array([[a, bb], [c, d]], float) + 0.5)
    st = StratifiedTable(tabs); lo, hi = st.oddsratio_pooled_confint(); return st.oddsratio_pooled, lo, hi
hs = hs.reset_index(drop=True); hs["band"] = hs.band.astype(int)
rows = []
for vn in list(FIXED) + list(LODO):
    for k in SETS:
        mset = FIXED[vn] if vn in FIXED else LODO[vn](k)
        if vn == "all 42": continue
        r = hs.motif.isin(mset).to_numpy(); r = pd.Series(r, index=hs.index); y = pd.Series(LAB[k], index=hs.index)
        o, lo, hi = mh(r, y)
        rows.append(dict(set=vn, heldout=k, n_motifs=len(mset), hits_kept=int(r.sum()), precision_kept=float(y[r].mean()), precision_removed=float(y[~r].mean()), mh_or=o, lo=lo, hi=hi))
HL = pd.DataFrame(rows); HL.to_csv(V4 / "v4_standard_hitlevel.tsv", sep="\t", index=False); pd.set_option("display.width", 220); print(HL.round(3).to_string(index=False))
hit_ok = {v: int((HL[HL.set == v].lo > 1).sum()) for v in HL.set.unique()}; print("hit-level sets with MH lower bound > 1 (of 6):", hit_ok)

# ---------------- (2) gene level: PPARgene verified human targets
G = pd.read_csv(f"{R}/benchmark/external_data/PPARgene_verified_targets.csv"); G = G[G.Species == "human"].copy(); G["sub"] = G.Subtype.str.upper(); G["SYM"] = G.Gene_symbol.str.upper()
sym = gh.SYMBOL.str.upper().to_numpy(); uniq = pd.Series(sym).map(pd.Series(sym).value_counts()).to_numpy() == 1
gold = {"human activation (primary)": set(G[G.Type == "Activation"].SYM), "human all": set(G.SYM), **{f"{s_} (all evidence)": set(G[G["sub"] == s_].SYM) for s_ in ["PPARA", "PPARD", "PPARG"]}}
for k, v in gold.items(): print(k, len(v), "in universe:", int(np.isin(sym[uniq], list(v)).sum()))
def score_for(motifs): return resid(raw_R(loci_table(1e-5, motifs=motifs), "noisy"))
def loci_weighted(wmap, gate=1e-5):
    hh = Hn[(Hn.p <= gate) & Hn.motif.isin(wmap.index[wmap > 0])].sort_values(["cidx", "c"]).copy()
    hh["locus"] = (((hh.cidx.diff().fillna(1) != 0) | (hh.c.diff().fillna(99) > 10))).cumsum().to_numpy(); hh["w"] = hh.motif.map(wmap)
    p = Pn[(Pn.d >= -3000) & (Pn.d <= 1000) & Pn.hit_id.isin(hh.hit_id)].merge(hh[["hit_id", "locus", "w"]], on="hit_id"); p["q"] = p.q * p.w
    lg = p.groupby(["gene_idx", "locus"]).q.max().reset_index(); lg["gp"] = gpos.reindex(lg.gene_idx).to_numpy(); return lg
wsoft = ((A.mean(axis=1).reindex(ALL) - 0.5) / 0.10).clip(0, 1).fillna(0)
X = {"all 42": score_for(None), "reliable-25": score_for(RELM_ALL), "PPAR-centric (9)": score_for(s2(RELM_ALL)), "no RXR half-sites (20)": score_for(s1(RELM_ALL)),
     "S0 remove RXR half-sites only (28)": score_for(s0), "DR1 (22)": score_for(dr1), "DR1-Q (20)": score_for(dr1q), "DR1-NR (5)": score_for(nr),
     "SOFT weights": resid(raw_R(loci_weighted(wsoft), "noisy"))}
def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
cv_lpl, cv_gc = cov.lpl.to_numpy(), cov.gc.to_numpy()
bins = pd.qcut(cv_lpl, 10, labels=False, duplicates="drop") * 10 + pd.qcut(cv_gc, 10, labels=False, duplicates="drop")
rg = np.random.default_rng(20261012)
def matched_auc(x, y, draws=50):
    pos = np.flatnonzero(y); out = []
    byb = {b: np.flatnonzero((bins == b) & ~y) for b in np.unique(bins)}
    for _ in range(draws):
        neg = np.concatenate([rg.choice(byb[bins[p]], 10, replace=True) for p in pos if len(byb[bins[p]])]); out.append(auc_of(np.r_[x[pos], x[neg]], np.r_[np.ones(len(pos), bool), np.zeros(len(neg), bool)]))
    return float(np.mean(out))
res = []; boot_pairs = {}
for gn, gs in gold.items():
    idx = np.flatnonzero(uniq); y = np.isin(sym[idx], list(gs)); pi, ni = np.flatnonzero(y), np.flatnonzero(~y)
    if y.sum() < 10: continue
    bs_idx = [(rg.choice(pi, len(pi)), rg.choice(ni, len(ni))) for _ in range(1000)]
    for vn, x in X.items():
        xi = x[idx]; a = auc_of(xi, y); ap = average_precision_score(y, xi); m = matched_auc(xi, y, 20)
        d = {"gold": gn, "set": vn, "n_pos": int(y.sum()), "auroc": a, "auprc": ap, "matched_auroc": m, "auprc_baseline": float(y.mean())}
        if vn != "all 42":
            x0 = X["all 42"][idx]; df = []
            for aa, bb in bs_idx:
                ii = np.r_[aa, bb]; yy = np.r_[np.ones(len(aa), bool), np.zeros(len(bb), bool)]; df.append(auc_of(xi[ii], yy) - auc_of(x0[ii], yy))
            d.update(d_auroc=a - auc_of(x0, y), d_lo=np.percentile(df, 2.5), d_hi=np.percentile(df, 97.5), d_auprc=ap - average_precision_score(y, x0), d_matched=m - matched_auc(x0, y, 20))
        res.append(d)
GL = pd.DataFrame(res); GL.to_csv(V4 / "v4_standard_gold.tsv", sep="\t", index=False); print(GL.round(3).to_string(index=False))

# ---------------- (3) window sweep, all 42, primary gold standard
idx = np.flatnonzero(uniq); y = np.isin(sym[idx], list(gold["human activation (primary)"]))
WS = []
for nm, (lo_, hi_) in {"-1000..+500": (-1000, 500), "-3000..+1000 (default)": (-3000, 1000), "-10000..+3000": (-10000, 3000)}.items():
    x = resid(raw_R(loci_table(1e-5, dlo=lo_, dhi=hi_), "noisy"))[idx]; WS.append(dict(window=nm, auroc=auc_of(x, y), auprc=average_precision_score(y, x), matched_auroc=matched_auc(x, y, 20)))
WS = pd.DataFrame(WS); WS.to_csv(V4 / "v4_standard_windows.tsv", sep="\t", index=False); print(WS.round(3).to_string(index=False))

# ---------------- rule
P_ = GL[(GL.gold == "human activation (primary)") & (GL.set != "all 42")]
for _, r in P_.iterrows():
    v = r.set; a_ = hit_ok.get(v, None)
    ra = (v == 'SOFT weights') or ((a_ is not None) and a_ >= 4);  # soft weights have no hit-level set: (a) not applicable
    rb = (r.d_lo > 0) and (r.d_auprc > 0) and (r.d_matched >= 0)
    print(f"{v}: (a) hit-level MH lower>1 in {a_}/6 -> {ra} | (b) AUROC diff {r.d_auroc:+.3f} ({r.d_lo:+.3f}), auPRC diff {r.d_auprc:+.4f}, matched diff {r.d_matched:+.3f} -> {rb} | ADOPT: {ra and rb}")
