"""Are high-scoring FIMO hits of the 42 motifs more often inside phastCons conserved elements than GC-matched control windows?  (no function labels)
Hits: p <= 1e-4 FIMO hits, non-repeat (UCSC RepeatMasker), grouped by motif class (heterodimer 3 / PPAR-only 11 / RXR-only 28) and by relative score (score / best score of the motif: <0.5, 0.5-0.6, ..., >=0.9).
Controls: one window per hit, same width, from a 1 kb window of the same GC class (2.5 % steps), N-free, non-repeat; for promoter hits (-3 kb/+1 kb of a TSS) the control comes from promoter windows (context-matched).
Conservation: UCSC phastCons 100-way conserved elements (hg38 phastConsElements100way, 5.2 % of the genome). A window is 'conserved' if >= 50 % of its bases are inside an element.
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats

rng = np.random.default_rng(30)
ROOT = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline"
FULL = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline_full_local"
CL = os.path.join(ROOT, "clustering")
FIMO = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar\FIMO_results_human42_unmasked"
INV = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar\Motifs_update_2026-09-30\motif_inventory_gates.csv"
CONS = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\R_code\benchmark\external_data\conservation\phastConsElements100way.txt.gz"
OUT = os.path.join(ROOT, "conservation"); os.makedirs(OUT, exist_ok=True)
CHRS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CI = {c: i for i, c in enumerate(CHRS)}
PER_CELL = 20000

def merge(s, e):
    o = np.argsort(s, kind="stable"); s, e = s[o], e[o]
    ce = np.maximum.accumulate(e); new = np.concatenate([[True], s[1:] > ce[:-1]]); idx = np.flatnonzero(new)
    return s[idx], np.maximum.reduceat(e, idx)

# ---- intervals ----
def read_bed(path, cols=(0, 1, 2), names=("c", "s", "e"), comp="infer"):
    d = pd.read_csv(path, sep="\t", header=None, usecols=list(cols), names=list(names), compression=comp)
    out = {}
    for c, g in d.groupby("c"):
        if c in CI: out[c] = merge(g.s.to_numpy(np.int64), g.e.to_numpy(np.int64))
    return out
reps = read_bed(os.path.join(CL, "hg38_repeats_merged.bed"))
cons = {}
d = pd.read_csv(CONS, sep="\t", header=None, usecols=[1, 2, 3], names=["c", "s", "e"], compression="gzip")
for c, g in d.groupby("c"):
    if c in CI: cons[c] = merge(g.s.to_numpy(np.int64), g.e.to_numpy(np.int64))
cons_cum = {c: (s, e, np.concatenate([[0], np.cumsum(e - s)])) for c, (s, e) in cons.items()}
print("conserved elements (merged):", sum(len(v[0]) for v in cons.values()), " covered Mb:", round(sum((v[1] - v[0]).sum() for v in cons.values()) / 1e6, 1), flush=True)

def hit_repeat(c, s, e):
    if c not in reps: return np.zeros(len(s), bool)
    st, en = reps[c]; k = np.searchsorted(st, e, "left") - 1
    return (k >= 0) & (en[np.clip(k, 0, None)] > s)

def cov_upto(c, x):
    st, en, cs = cons_cum[c]; k = np.searchsorted(st, x, "right") - 1
    kk = np.clip(k, 0, None); inside = np.clip(x - st[kk], 0, en[kk] - st[kk])
    return np.where(k >= 0, cs[kk] + inside, 0)

def overlap_frac(c, s, e):
    if c not in cons_cum: return np.zeros(len(s))
    return (cov_upto(c, e) - cov_upto(c, s)) / (e - s)

# ---- GC windows and promoter windows ----
gc = pd.read_csv(os.path.join(CL, "gc_bins_1kb.tsv"), sep="\t")
gc = gc[gc.chrom.isin(CHRS)]
gc["gcf"] = gc.gc / (gc.length - gc.n).clip(lower=1); gc["cls"] = np.floor(gc.gcf * 100 / 2.5).astype(int)
gc["ok"] = (gc.n == 0) & (gc.length == 1000)
cls_map = {}; maxbin = {}
for c, g in gc.groupby("chrom"):
    arr = np.full(g.bin.max() + 1, -1, np.int32); arr[g.bin.to_numpy()] = np.where(g.ok.to_numpy(), g.cls.to_numpy(), -1); cls_map[c] = arr
tss = pd.read_csv(os.path.join(FULL, "tss_universe.tsv"), sep="\t")
prom_tss = {}
for (c, st), g in tss.groupby(["chr", "strand"]):
    if c in CI: prom_tss[(c, st)] = np.sort(g.tss.to_numpy(np.int64))
def in_promoter(c, mid):
    res = np.zeros(len(mid), bool)
    for st, lo, hi in (("+", -1000, 3000), ("-", -3000, 1000)):   # TSS must lie in [mid+lo, mid+hi]
        t = prom_tss.get((c, st))
        if t is None: continue
        res |= (np.searchsorted(t, mid + hi, "right") - np.searchsorted(t, mid + lo, "left")) > 0
    return res
prom_bins = {}
for c in CHRS:
    nb = len(cls_map[c]); m = np.zeros(nb, bool)
    for (cc, st), t in prom_tss.items():
        if cc != c: continue
        lo = (t - (3000 if st == "+" else 1000)) // 1000; hi = (t + (1000 if st == "+" else 3000)) // 1000
        for a, b in zip(lo, hi):
            m[max(a, 0):min(b + 1, nb)] = True
    prom_bins[c] = m
# eligible window lists per GC class: genome-wide and promoter
elig = {"all": {}, "prom": {}}
for c in CHRS:
    cm = cls_map[c]; ok = cm >= 0
    for key, mask in (("all", ok), ("prom", ok & prom_bins[c])):
        for k in np.unique(cm[mask]):
            elig[key].setdefault(int(k), []).append((c, np.flatnonzero(mask & (cm == k))))
for key in elig:
    for k, lst in elig[key].items():
        elig[key][k] = (np.concatenate([np.full(len(b), CI[c], np.int16) for c, b in lst]), np.concatenate([b for c, b in lst]))
print("GC classes with windows:", len(elig["all"]), " promoter:", len(elig["prom"]), flush=True)

def controls(df, key):
    """df: chrom,start,stop,width,cls. one GC-class-matched non-repeat control window per row."""
    n = len(df); cc = np.full(n, -1, np.int16); cs = np.zeros(n, np.int64)
    todo = np.arange(n)
    for _ in range(12):
        if len(todo) == 0: break
        new_c = np.zeros(len(todo), np.int16); new_s = np.zeros(len(todo), np.int64); bad = np.zeros(len(todo), bool)
        for k in np.unique(df.cls.to_numpy()[todo]):
            sel = np.flatnonzero(df.cls.to_numpy()[todo] == k)
            if k not in elig[key]: bad[sel] = True; continue
            ch, bn = elig[key][k]; pick = rng.integers(0, len(ch), len(sel))
            w = df.width.to_numpy()[todo][sel]; off = (rng.random(len(sel)) * (1000 - w)).astype(np.int64)
            new_c[sel] = ch[pick]; new_s[sel] = bn[pick].astype(np.int64) * 1000 + off
        w = df.width.to_numpy()[todo]; rep = np.zeros(len(todo), bool)
        for ci in np.unique(new_c):
            m = new_c == ci; rep[m] = hit_repeat(CHRS[ci], new_s[m], new_s[m] + w[m])
        good = ~rep & ~bad; cc[todo[good]] = new_c[good]; cs[todo[good]] = new_s[good]; todo = todo[~good]
    ok = cc >= 0
    return cc, cs, ok

# ---- hits ----
inv = pd.read_csv(INV, encoding="utf-8-sig"); catmap = dict(zip(inv.motif_id, inv.category))
grp = {"PPAR:RXR heterodimer": "heterodimer 3", "PPAR": "PPAR-only 11", "RXR": "RXR-only 28"}
bks = [(-1, 0.5, "<0.5"), (0.5, 0.6, "0.5-0.6"), (0.6, 0.7, "0.6-0.7"), (0.7, 0.8, "0.7-0.8"), (0.8, 0.9, "0.8-0.9"), (0.9, 9, ">=0.9")]
samples = []
for b in sorted(glob.glob(os.path.join(FIMO, "*", "fimo.tsv"))):
    d = pd.read_csv(b, sep="\t", usecols=["motif_id", "sequence_name", "start", "stop", "score"], dtype={"sequence_name": str}, comment="#")
    d = d[d.sequence_name.isin([str(i) for i in range(1, 23)] + ["X"])].copy()
    d["chrom"] = "chr" + d.sequence_name; d["rel"] = d.score / d.groupby("motif_id").score.transform("max")
    d["group"] = d.motif_id.map(catmap).map(grp); d["start"] = d.start.astype(np.int64) - 1; d["stop"] = d.stop.astype(np.int64)
    for lo, hi, name in bks:
        for g in grp.values():
            s = d[(d.group == g) & (d.rel >= lo) & (d.rel < hi)]
            if len(s): samples.append(s.sample(min(len(s), PER_CELL), random_state=int(rng.integers(1 << 30))).assign(bucket=name)[["chrom", "start", "stop", "group", "bucket", "rel", "motif_id"]])
    print(os.path.basename(os.path.dirname(b)), "sampled", flush=True)
H = pd.concat(samples, ignore_index=True)
rep = np.zeros(len(H), bool)
for c, ix in H.groupby("chrom").indices.items(): rep[ix] = hit_repeat(c, H.start.to_numpy()[ix], H.stop.to_numpy()[ix])
H = H[~rep].copy()
H = H.groupby(["group", "bucket"], group_keys=False).apply(lambda x: x.sample(min(len(x), PER_CELL), random_state=1)).reset_index(drop=True)
H["width"] = H.stop - H.start; H["mid"] = (H.start + H.stop) // 2
H["prom"] = False
for c, ix in H.groupby("chrom").indices.items(): H.loc[H.index[ix], "prom"] = in_promoter(c, H.mid.to_numpy()[ix])
cl = np.full(len(H), -1)
for c, ix in H.groupby("chrom").indices.items():
    b = H.mid.to_numpy()[ix] // 1000; arr = cls_map[c]; cl[ix] = np.where(b < len(arr), arr[np.clip(b, 0, len(arr) - 1)], -1)
H["cls"] = cl; H = H[H.cls >= 0].reset_index(drop=True)
print("hits kept:", len(H), " in promoters:", int(H.prom.sum()), flush=True)
H["hit_ov"] = 0.0
for c, ix in H.groupby("chrom").indices.items(): H.loc[H.index[ix], "hit_ov"] = overlap_frac(c, H.start.to_numpy()[ix], H.stop.to_numpy()[ix])
res = []
for analysis, mask, key in (("genome-wide (non-repeat)", np.ones(len(H), bool), "all"), ("promoter hits vs promoter controls", H.prom.to_numpy(), "prom")):
    D = H[mask].reset_index(drop=True); cc, cs, ok = controls(D, key); D = D[ok].reset_index(drop=True); cc, cs = cc[ok], cs[ok]
    ov = np.zeros(len(D))
    for ci in np.unique(cc):
        m = cc == ci; ov[m] = overlap_frac(CHRS[ci], cs[m], cs[m] + D.width.to_numpy()[m])
    D["ctrl_ov"] = ov
    for (g, bk), x in D.groupby(["group", "bucket"]):
        h = (x.hit_ov >= 0.5).to_numpy(); c = (x.ctrl_ov >= 0.5).to_numpy()
        n = len(x); ph, pc = h.mean(), c.mean()
        table = [[h.sum(), n - h.sum()], [c.sum(), n - c.sum()]]
        orr, p = stats.fisher_exact(table) if n < 100000 else (np.nan, np.nan)
        bs = [(rng.choice(h, n).mean() - rng.choice(c, n).mean()) for _ in range(200)]
        res.append(dict(analysis=analysis, group=g, bucket=bk, n=n, frac_conserved_hits=ph, frac_conserved_controls=pc, diff=ph - pc, diff_lo=np.percentile(bs, 2.5), diff_hi=np.percentile(bs, 97.5),
                        fold=ph / pc if pc > 0 else np.nan, odds_ratio=orr, fisher_p=p, mean_overlap_hits=x.hit_ov.mean(), mean_overlap_controls=x.ctrl_ov.mean()))
    print(analysis, "done:", len(D), "pairs", flush=True)
R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "conservation_of_hits_vs_gc_matched_controls.csv"), index=False)
pd.set_option("display.width", 250)
print(R[["analysis", "group", "bucket", "n", "frac_conserved_hits", "frac_conserved_controls", "fold", "diff", "diff_lo", "diff_hi"]].round(4).to_string(index=False))
