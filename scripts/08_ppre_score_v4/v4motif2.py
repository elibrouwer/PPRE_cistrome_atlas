"""Dropping low-quality and most RXR motifs (DESIGN_V4.md Addendum 10)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))
MT = pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_gc_band_auc_by_set.tsv", sep="\t")
SM = pd.read_csv(V4 / "motif_summary_42.csv").set_index("motif")
hit_motifs = set(Hn.motif); ALL = sorted(hit_motifs)
A = MT.groupby(["motif", "set"]).auc.mean().unstack(); O = MT[MT.split == "odd"].set_index(["motif", "set"])
def select_excluding(k):
    cols = [c for c in A.columns if c != k]; sub = A[cols]; nh = O.n_hits.unstack()[cols]; npz = O.n_pos.unstack()[cols]
    ok = (sub.median(axis=1) >= 0.55) & ((sub > 0.5).sum(axis=1) >= 4) & (nh.median(axis=1) >= 2000) & (npz.median(axis=1) >= 50)
    return set(ok[ok].index) & hit_motifs
def s2(mset): return {m for m in mset if SM.loc[m, "type"] in ("PPAR", "heterodimer")}
def s1(mset): return {m for m in mset if not (SM.loc[m, "type"] == "RXR" and SM.loc[m, "length_bp"] <= 14)}
VAR = {"S2 PPAR-centric": s2, "S1 no RXR half-sites": s1, "reliable-25 (ref)": lambda m: m}
def score_for(motifs): return resid(raw_R(loci_table(1e-5, motifs=motifs), "noisy"))
x_all = score_for(None); rs = np.random.default_rng(20261010); rows = []
for vn, fn in VAR.items():
    for k in SETS:
        sel = fn(select_excluding(k)); xs = score_for(sel); o_all, o_sel = adj(CHIPY[k], x_all), adj(CHIPY[k], xs)
        rnd = [adj(CHIPY[k], score_for(set(rs.choice(ALL, len(sel), replace=False)))) for _ in range(30)]
        rows.append(dict(variant=vn, heldout=k, n_selected=len(sel), or_all42=o_all, or_selected=o_sel, random_p90=np.percentile(rnd, 90), beats_all=o_sel >= o_all, beats_p90=o_sel > np.percentile(rnd, 90)))
    say(f"{vn} LODO done")
T = pd.DataFrame(rows); T.to_csv(V4 / "v4motif2_lodo_chip.tsv", sep="\t", index=False); print(T.round(3).to_string(index=False))
full = select_excluding.__wrapped__ if hasattr(select_excluding, "__wrapped__") else None
RELM_ALL = set(pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_reliability.tsv", sep="\t").query("reliable").motif) & hit_motifs
lit0 = auc_of(x_all, ylit); out = []
for vn, fn in VAR.items():
    mset = fn(RELM_ALL); xv = score_for(mset); hr = []
    for nm, idx, y in CON:
        pi, ni = np.flatnonzero(y), np.flatnonzero(~y); a0, a1 = x_all[idx], xv[idx]; df = []
        for _ in range(500):
            a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
            df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
        hr.append(dict(variant=vn, n_motifs=len(mset), contrast=nm, auc_all42=auc_of(a0, y), auc_variant=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
    H_ = pd.DataFrame(hr); out.append(H_); lit1 = auc_of(xv, ylit)
    ca, cb = int(T[T.variant == vn].beats_all.sum()), int(T[T.variant == vn].beats_p90.sum()); ra = ca >= 4 and cb >= 3; rb = int((H_.lo > 0).sum()) >= 3; rc = lit1 >= lit0 - 0.01
    print(f"{vn}: motifs {len(mset)} | (a) beats all-42 {ca}/6, beats random p90 {cb}/6 -> {ra} | (b) hepatocyte lower>0 {int((H_.lo > 0).sum())}/4 -> {rb} | (c) literature {lit1:.3f} vs {lit0:.3f} -> {rc} | ADOPT: {ra and rb and rc}")
R_ = pd.concat(out); R_.to_csv(V4 / "v4motif2_hepatocyte.tsv", sep="\t", index=False); print(R_.round(3).to_string(index=False))
