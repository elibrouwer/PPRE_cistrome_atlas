"""Repeat-embedded hits and motif AUC (DESIGN_V4.md Addendum 11)."""
from pathlib import Path
import pyarrow.parquet as pq
src = Path("v4_windows.py").read_text(encoding="utf8")
exec(compile(src[:src.index("# ---- liver ATAC")], "w", "exec"))
SM = pd.read_csv(V4 / "motif_summary_42.csv").set_index("motif")
chrom_ix = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
rm = pd.read_csv(Path("C:/Users/brouw/rmsk_hg38.bed"), sep="\t", header=None, usecols=[0, 1, 2]); rm.columns = ["c", "s", "e"]
rm = rm[rm.c.isin(chrom_ix)].copy(); rm["ci"] = rm.c.map(chrom_ix); rm = rm.sort_values(["ci", "s"]).reset_index(drop=True)
rm["pm"] = rm.groupby("ci").e.cummax().groupby(rm.ci).shift().fillna(-1)
g = ((rm.ci.diff().fillna(1) != 0) | (rm.s > rm.pm)).cumsum(); rm = rm.groupby(g).agg(ci=("ci", "first"), s=("s", "min"), e=("e", "max")).reset_index(drop=True)
RMS = IntervalSet(rm.ci.to_numpy(), rm.s.to_numpy(), rm.e.to_numpy()); say(f"repeat intervals {len(rm):,}")
H = pd.read_parquet(CACHE / "human42" / "hits.parquet", columns=["motif_id", "cidx", "start0", "end", "pvalue", "flank_gc"])
H["motif_id"] = H.motif_id.astype(str); H = H[H.cidx % 2 == 1].reset_index(drop=True)
H["rep"] = RMS.overlaps(H.cidx.to_numpy(), H.start0.to_numpy(), H.end.to_numpy()); H["s"] = -np.log10(H.pvalue.clip(lower=1e-12)); say(f"hits (odd chromosomes) {len(H):,}, in repeat {H.rep.mean():.3f}")
for s_ in SETS: H[s_] = PK[s_].overlaps(H.cidx.to_numpy(), H.start0.to_numpy(), H.end.to_numpy())
def band_auc(d, y):
    if y.sum() < 20 or (~y).sum() < 20: return np.nan
    b = pd.qcut(d.flank_gc, 10, labels=False, duplicates="drop"); num = den = 0
    for k in np.unique(b):
        m = (b == k).to_numpy(); yy = y[m]; n1 = yy.sum(); n0 = (~yy).sum()
        if n1 < 5 or n0 < 5: continue
        r = stats.rankdata(d.s.to_numpy()[m]); a = (r[yy].sum() - n1 * (n1 + 1) / 2) / (n1 * n0); w = n1 * n0; num += a * w; den += w
    return num / den if den else np.nan
rows = []
for m, d in H.groupby("motif_id"):
    r = dict(motif=m, repeat_share=d.rep.mean(), top10_repeat_share=d[d.s >= d.s.quantile(0.9)].rep.mean(), reliable=bool(SM.loc[m, "reliable_25"]))
    for lab, sub in [("all", d), ("nonrepeat", d[~d.rep]), ("repeat", d[d.rep])]:
        r[f"auc_{lab}"] = np.nanmedian([band_auc(sub, sub[s_].to_numpy()) for s_ in SETS])
    rows.append(r)
T = pd.DataFrame(rows).set_index("motif"); T.to_csv(V4 / "v4repeats_motifs.tsv", sep="\t"); pd.set_option("display.width", 200)
low = ["MA0066.1", "MA0066.2", "PPARG_M02901_3.10", "RXRB_M10655_3.10"]
print(T.loc[low].round(3).to_string()); med = T[T.reliable].repeat_share.median(); print("median repeat share of the 25 reliable motifs:", round(med, 3), "| their median auc all / nonrepeat:", T[T.reliable].auc_all.median().round(3), T[T.reliable].auc_nonrepeat.median().round(3))
ok = [(T.loc[m, "repeat_share"] > med) and (T.loc[m, "auc_nonrepeat"] - T.loc[m, "auc_all"] >= 0.03) for m in low]
print(dict(zip(low, ok)), "-> repeat explanation supported:", sum(ok) >= 3)
