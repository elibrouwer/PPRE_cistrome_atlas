"""Leave-one-dataset-out motif prioritisation (DESIGN_V4.md Addendum 5)."""
from pathlib import Path
src = Path("v4robust.py").read_text(encoding="utf8")
exec(compile(src[:src.index("rows = []\nfor nm, x in SC.items()")], "r", "exec"))   # hits, loci_table, raw_R, resid, CON, ylit, CHIPY, adj, auc_of, RELM
MT = pd.read_csv(f"{R}/results/PPRE_score_v3/side_reliable_motifs/tables/motif_gc_band_auc_by_set.tsv", sep="\t")
hit_motifs = set(Hn.motif); ALL = sorted(hit_motifs); print("motifs in hits:", len(ALL), "| in AUC table:", len(set(MT.motif) & hit_motifs))
A = MT.groupby(["motif", "set"]).auc.mean().unstack()            # mean of odd and even
O = MT[MT.split == "odd"].set_index(["motif", "set"])
def select_excluding(k):
    cols = [c for c in A.columns if c != k]; sub = A[cols]
    nh = O.n_hits.unstack()[cols]; npz = O.n_pos.unstack()[cols]
    ok = (sub.median(axis=1) >= 0.55) & ((sub > 0.5).sum(axis=1) >= 4) & (nh.median(axis=1) >= 2000) & (npz.median(axis=1) >= 50)
    return set(ok[ok].index) & hit_motifs
def score_for(motifs): return resid(raw_R(loci_table(1e-5, motifs=motifs), "noisy"))
x_all = score_for(None)
rs = np.random.default_rng(20261009)
rows = []
for k in SETS:
    sel = select_excluding(k); xs = score_for(sel)
    o_all, o_sel = adj(CHIPY[k], x_all), adj(CHIPY[k], xs)
    rnd = [adj(CHIPY[k], score_for(set(rs.choice(ALL, len(sel), replace=False)))) for _ in range(30)]
    rows.append(dict(heldout=k, n_selected=len(sel), or_all42=o_all, or_selected=o_sel, random_mean=np.mean(rnd), random_p90=np.percentile(rnd, 90), beats_all=o_sel >= o_all, beats_p90=o_sel > np.percentile(rnd, 90)))
    say(f"held-out {k}: {len(sel)} motifs, OR {o_sel:.3f} vs all-42 {o_all:.3f}, random p90 {np.percentile(rnd, 90):.3f}")
T = pd.DataFrame(rows); T.to_csv(V4 / "v4motif_lodo_chip.tsv", sep="\t", index=False); print(T.round(3).to_string(index=False))
ca, cb = int(T.beats_all.sum()), int(T.beats_p90.sum()); rule_a = ca >= 4 and cb >= 3
# (b) hepatocyte with the 25-motif set chosen on all 6 sets
x25 = score_for(RELM); hr = []
for nm, idx, y in CON:
    pi, ni = np.flatnonzero(y), np.flatnonzero(~y); a0, a1 = x_all[idx], x25[idx]; df = []
    for _ in range(500):
        a = rng.choice(pi, len(pi)); b = rng.choice(ni, len(ni)); ii = np.r_[a, b]; yy = np.r_[np.ones(len(a), bool), np.zeros(len(b), bool)]
        df.append(auc_of(a1[ii], yy) - auc_of(a0[ii], yy))
    hr.append(dict(contrast=nm, auc_all42=auc_of(a0, y), auc_25=auc_of(a1, y), diff=auc_of(a1, y) - auc_of(a0, y), lo=np.percentile(df, 2.5), hi=np.percentile(df, 97.5)))
H_ = pd.DataFrame(hr); H_.to_csv(V4 / "v4motif_hepatocyte.tsv", sep="\t", index=False); print(H_.round(3).to_string(index=False))
rule_b = int((H_.lo > 0).sum()) >= 3
lit0, lit1 = auc_of(x_all, ylit), auc_of(x25, ylit); rule_c = lit1 >= lit0 - 0.01
print(f"(a) held-out ChIP: beats all-42 in {ca}/6 (need 4), beats random p90 in {cb}/6 (need 3) -> {rule_a}")
print(f"(b) hepatocyte lower bound > 0 in {int((H_.lo > 0).sum())}/4 (need 3) -> {rule_b}")
print(f"(c) literature AUROC all-42 {lit0:.3f}, 25-motif {lit1:.3f} (need >= -0.01) -> {rule_c}")
print("MOTIF PRIORITISATION WORKS:", rule_a and rule_b and rule_c)
