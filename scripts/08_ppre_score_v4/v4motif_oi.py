"""Outcome-independent motif sets and soft motif weights (DESIGN_V4.md Addendum 14)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
INV = pd.read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/motif_inventory_human_only.csv", encoding="utf-8-sig").set_index("motif_id")
MT = pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_gc_band_auc_by_set.tsv", sep="\t")
A = MT.groupby(["motif", "set"]).auc.mean().unstack()
hit_motifs = set(Hn.motif); ALL = sorted(hit_motifs)
dr1 = set(INV.index[INV.architecture == "DR1"]) & hit_motifs
dr1q = {m for m in dr1 if str(INV.loc[m, "sparse_flag"]) != "SPARSE"}
sub = INV.loc[sorted(dr1q)].assign(wid=INV.loc[sorted(dr1q), "width"]).sort_values(["ic_bits", "wid"], ascending=False)
nr = set(sub.groupby("cluster_r098").head(1).index)
SETS_OI = {"DR1": dr1, "DR1-Q": dr1q, "DR1-NR (primary)": nr}
for k, v in SETS_OI.items(): print(k, len(v), sorted(v))
def score_for(motifs): return resid(raw_R(loci_table(1e-5, motifs=motifs), "noisy"))
def loci_weighted(wmap, gate=1e-5):
    hs = Hn[(Hn.p <= gate) & Hn.motif.isin(wmap.index[wmap > 0])].sort_values(["cidx", "c"]).copy()
    hs["locus"] = (((hs.cidx.diff().fillna(1) != 0) | (hs.c.diff().fillna(99) > 10))).cumsum().to_numpy(); hs["w"] = hs.motif.map(wmap)
    p = Pn[(Pn.d >= -3000) & (Pn.d <= 1000) & Pn.hit_id.isin(hs.hit_id)].merge(hs[["hit_id", "locus", "w"]], on="hit_id"); p["q"] = p.q * p.w
    lg = p.groupby(["gene_idx", "locus"]).q.max().reset_index(); lg["gp"] = gpos.reindex(lg.gene_idx).to_numpy(); return lg
def soft_score(sets_for_auc):
    auc = A[sets_for_auc].mean(axis=1).reindex(ALL); w = ((auc - 0.5) / 0.10).clip(0, 1).fillna(0); return resid(raw_R(loci_weighted(w), "noisy")), w
x_all = score_for(None); rs = np.random.default_rng(20261011); rows = []
# (a) outcome-independent sets: random control sampled once per variant
for vn, mset in SETS_OI.items():
    xs = score_for(mset); rnd = [score_for(set(rs.choice(ALL, len(mset), replace=False))) for _ in range(30)]
    for k in SETS:
        o_all, o_sel = adj(CHIPY[k], x_all), adj(CHIPY[k], xs); p90 = np.percentile([adj(CHIPY[k], r_) for r_ in rnd], 90)
        rows.append(dict(variant=vn, heldout=k, n=len(mset), or_all42=o_all, or_selected=o_sel, random_p90=p90, beats_all=o_sel >= o_all, beats_p90=o_sel > p90))
    say(f"{vn} done")
# SOFT with leave-one-dataset-out weights, random control = permuted weights
for k in SETS:
    xs, w = soft_score([c for c in SETS if c != k]); o_all, o_sel = adj(CHIPY[k], x_all), adj(CHIPY[k], xs)
    rnd = []
    for _ in range(30):
        wp = pd.Series(rs.permutation(w.values), index=w.index); rnd.append(adj(CHIPY[k], resid(raw_R(loci_weighted(wp), "noisy"))))
    rows.append(dict(variant="SOFT", heldout=k, n=int((w > 0).sum()), or_all42=o_all, or_selected=o_sel, random_p90=np.percentile(rnd, 90), beats_all=o_sel >= o_all, beats_p90=o_sel > np.percentile(rnd, 90)))
say("SOFT LODO done")
T = pd.DataFrame(rows); T.to_csv(V4 / "v4motif_oi_chip.tsv", sep="\t", index=False); print(T.round(3).to_string(index=False))
lit0 = auc_of(x_all, ylit); out = []
variants = {**{k: score_for(v) for k, v in SETS_OI.items()}, "SOFT": soft_score(SETS)[0]}
for vn, xv in variants.items():
    hr = []
    for nm, idx, y in CON:
        pi, ni = np.flatnonzero(y), np.flatnonzero(~y); a0, a1 = x_all[idx], xv[idx]; df = []
        for _ in range(500):
            a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
            df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
        hr.append(dict(variant=vn, contrast=nm, auc_all42=auc_of(a0, y), auc_variant=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
    H_ = pd.DataFrame(hr); out.append(H_); lit1 = auc_of(xv, ylit); Tv = T[T.variant == vn]
    ca, cb = int(Tv.beats_all.sum()), int(Tv.beats_p90.sum()); ra = ca >= 4 and cb >= 3; rb = int((H_.lo > 0).sum()) >= 3; rc = lit1 >= lit0 - 0.01
    print(f"{vn}: (a) beats all-42 {ca}/6, random p90 {cb}/6 -> {ra} | (b) hepatocyte lower>0 {int((H_.lo > 0).sum())}/4 -> {rb} | (c) literature {lit1:.3f} vs {lit0:.3f} -> {rc} | ADOPT: {ra and rb and rc}")
R_ = pd.concat(out); R_.to_csv(V4 / "v4motif_oi_hepatocyte.tsv", sep="\t", index=False); print(R_.round(3).to_string(index=False))
