"""PPAR/RXR-specificity of the hit-level binding signal: real 42 motifs vs column-permuted and reversed copies (DESIGN_V4.md Addendum 17).
Run from results/PPRE_score_v4 with env PPRE_CACHE_DIR, HPEAKS and CTRL42_OUT (folder with real/, perm01..05/, rev/ FIMO outputs scanned on the promoter union)."""
import os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
V4 = Path(__file__).resolve().parent if "__file__" in globals() else Path(".")
V4 = V4.parent if V4.name == "ctrl42" else V4
HERE = V4.parent / "PPRE_score_v3"
CACHE = Path(os.environ["PPRE_CACHE_DIR"]); CTRL = Path(os.environ["CTRL42_OUT"])
src = (HERE / "v3_human.py").read_text(encoding="utf8")
_own = globals().get("__file__"); __file__ = str(HERE / "v3_human.py")
exec(compile(src[:src.index("# ---- hit level")], "head", "exec"))                  # PK, SETS, label, say, ...
__file__ = _own
sys.path.insert(0, str(V4.parent / "PPARA_liver_ChIP_validation" / "scripts"))
from annot import load_gc_bins, LocalGC          # noqa
from common import STD_CHR                       # noqa
lg = LocalGC(load_gc_bins("human"))
OUT = V4 / "ctrl42"; OUT.mkdir(exist_ok=True)
ENS = {str(i): i - 1 for i in range(1, 23)}; ENS.update({"X": 22, "Y": 23})

def read_hits(setdir):
    parts = []
    for ch in pd.read_csv(setdir / "fimo.tsv", sep="\t", comment="#", chunksize=3_000_000, dtype={"sequence_name": str, "motif_id": str},
                          usecols=["motif_id", "sequence_name", "start", "stop", "p-value"]):
        sp = ch.sequence_name.str.split("|", expand=True)
        ch["cidx"] = sp[0].astype(int); w0 = sp[1].astype(np.int64)
        ch["start0"] = (w0 + ch.start - 1).astype(np.int64); ch["end"] = (w0 + ch.stop).astype(np.int64)
        parts.append(ch[["motif_id", "cidx", "start0", "end", "p-value"]].rename(columns={"p-value": "pvalue"}))
    H = pd.concat(parts, ignore_index=True)
    H["motif"] = H.motif_id.str.replace(r"__(P\d+|REV)$", "", regex=True)
    cen = (H.start0.to_numpy() + H.end.to_numpy()) // 2
    H["flank_gc"] = lg.query(STD_CHR, H.cidx.to_numpy(), cen, half=250)
    H["s"] = (-np.log10(H.pvalue.clip(lower=1e-12))).clip(4, 12)
    return H

def strat_auc(y, sc, gcv, nb=10):
    q = pd.qcut(gcv, nb, labels=False, duplicates="drop"); num = den = 0.0
    for b in np.unique(q[~pd.isna(q)]):
        m = q == b; yy = y[m]
        if yy.min() == yy.max(): continue
        w = yy.sum() * (~yy).sum(); num += w * roc_auc_score(yy, sc[m]); den += w
    return num / den if den > 0 else np.nan

rng = np.random.default_rng(20261010); rows = []; counts = []
for setname in ["real", "perm01", "perm02", "perm03", "perm04", "perm05", "rev"]:
    t0 = time.time(); H = read_hits(CTRL / setname); say(f"{setname}: {len(H):,} hits, {H.motif.nunique()} motifs")
    LAB = {s_: label(H, PK[s_]).astype(bool) for s_ in SETS}
    for mot, g in H.groupby("motif"):
        counts.append(dict(set=setname, motif=mot, n_hits=len(g)))
        idx = g.index.to_numpy(); gcv = g.flank_gc.to_numpy(); sc = g.s.to_numpy()
        for s_ in SETS:
            y = LAB[s_][idx]
            if y.sum() < 10 or (~y).sum() < 10: rows.append(dict(set=setname, motif=mot, chip=s_, n_pos=int(y.sum()), auc=np.nan)); continue
            pos = np.flatnonzero(y); neg = np.flatnonzero(~y)
            if len(neg) > 300000: neg = rng.choice(neg, 300000, replace=False)
            ii = np.r_[pos, neg]; yy = np.r_[np.ones(len(pos), bool), np.zeros(len(neg), bool)]
            rows.append(dict(set=setname, motif=mot, chip=s_, n_pos=int(len(pos)), auc=strat_auc(yy, sc[ii], gcv[ii])))
    say(f"  {setname} done in {time.time()-t0:.0f}s")
A = pd.DataFrame(rows); A.to_csv(OUT / "ctrl42_auc_by_chip.tsv", sep="\t", index=False)
Cn = pd.DataFrame(counts); Cn.to_csv(OUT / "ctrl42_hit_counts.tsv", sep="\t", index=False)
# pipeline check: real rescan counts vs genome-wide scan restricted to the promoter union
H42 = pd.read_parquet(CACHE / "human42" / "hits.parquet"); gw = H42.groupby(H42.motif_id.astype(str)).size()
chk = Cn[Cn.set == "real"].set_index("motif").n_hits.to_frame("rescan").join(gw.rename("genome_scan")); chk["ratio"] = chk.rescan / chk.genome_scan
say(f"pipeline check: rescan/genome-scan promoter hit ratio median {chk.ratio.median():.4f}, range {chk.ratio.min():.4f}-{chk.ratio.max():.4f}; within 1% for {(chk.ratio.sub(1).abs() <= 0.01).sum()} of {len(chk)} motifs")
M = A.groupby(["set", "motif"]).auc.median().unstack(0)
perm = M[[f"perm0{k}" for k in range(1, 6)]].mean(axis=1); M["perm_mean"] = perm
M["real_minus_perm"] = M.real - perm; M["real_minus_rev"] = M.real - M.rev; M["real_above_all_6"] = M.real > M[[f"perm0{k}" for k in range(1, 6)] + ["rev"]].max(axis=1)
M.to_csv(OUT / "ctrl42_motif_summary.tsv", sep="\t")
a = stats.wilcoxon(M.real_minus_perm.dropna(), alternative="greater"); b = (M.real_minus_perm > 0).mean()
c = stats.wilcoxon(M.real_minus_rev.dropna(), alternative="greater"); d_ = (M.real_minus_rev > 0).mean()
say(f"REAL median AUC {M.real.median():.3f}; permuted mean-of-means {perm.median():.3f}; reversed {M.rev.median():.3f}")
say(f"(a) real > permuted: signed-rank p = {a.pvalue:.3g} -> {a.pvalue < 0.05} | (b) share of motifs real > permuted mean {b:.2f} -> {b >= 0.6} | (c) real > reversed: p = {c.pvalue:.3g}, share {d_:.2f} -> {c.pvalue < 0.05}")
say(f"motifs with real above all 6 controls: {int(M.real_above_all_6.sum())} of {len(M)}")
spec = (a.pvalue < 0.05) and (b >= 0.6)
print("RESULT:", "specific against permutation and reversal" if (spec and c.pvalue < 0.05) else ("specific against permutation only" if spec else "NOT specific (composition effect)"))
