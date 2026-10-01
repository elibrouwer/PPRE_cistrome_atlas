"""Cluster propensity (CP) of the FIMO motif sites, after Wang, Wang & Zang, NAR 2025 (53:gkaf015), for the 42 motifs.

For every motif:
  * distance from each site to the next non-overlapping site downstream (same chromosome); log10, clipped at 1 bp
  * control sites: the same number of sites with the same width, placed uniformly at random
      A  "paper"  control: anywhere on chr1-22,X outside the assembly N gaps
      B  "repeat" control: keeps the repeat / non-repeat proportion of the real sites (sites inside RepeatMasker elements are
         placed uniformly inside repeats, the others uniformly outside), so repeat-driven clustering is not counted as clustering
      C  "GC-matched" control: keeps the GC distribution of the real sites; the genome is cut into 1 kb windows, windows are
         grouped into 40 GC classes, and control sites are drawn from windows of the same GC class as the real sites,
         so clustering that only reflects GC-rich regions (CpG islands, GC-rich isochores) is not counted as clustering
  * CP = signed two-sample Kolmogorov-Smirnov statistic between the observed and the control distance distributions
         (positive: observed sites are closer to their neighbours than control sites), averaged over n_iter random controls
The ENCODE hg38 blacklist v2 (hg38-blacklist.v2.bed.gz) is excluded by default: sites inside it are dropped and control sites
never fall in it (--no-blacklist switches this off). Not done here (compared with the paper): 100 iterations (default 10).

Usage:  python 11b_motif_cluster_propensity.py [--iter 10] [--tier 1e-4] [--motifs ID1,ID2]
"""
import argparse, glob, os, time
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

ROOT = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar"
FIMO = os.path.join(ROOT, "FIMO_results_human42_unmasked")
INV = os.path.join(ROOT, "Motifs_update_2026-09-30", "motif_inventory_human_only.csv")
CL = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline\clustering"

ap = argparse.ArgumentParser()
ap.add_argument("--iter", type=int, default=10)
ap.add_argument("--tier", type=float, default=1e-4)
ap.add_argument("--motifs", default="")
ap.add_argument("--no-blacklist", action="store_true")
args = ap.parse_args()

chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
cidx = {c: i for i, c in enumerate(chroms)}
SHIFT = np.int64(1) << 32


def read_bed(path):
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
    return {c: (g["s"].to_numpy(np.int64), g["e"].to_numpy(np.int64)) for c, g in df.groupby("c")}


sizes = dict(pd.read_csv(os.path.join(CL, "hg38_chrom_sizes.tsv"), sep="\t", header=None).values)
gaps, reps = read_bed(os.path.join(CL, "hg38_N_gaps.bed")), read_bed(os.path.join(CL, "hg38_repeats_merged.bed"))
bl = {} if args.no_blacklist else read_bed(os.path.join(CL, "hg38-blacklist.v2.bed.gz"))
EMPTY = (np.array([], np.int64), np.array([], np.int64))


def merge(s, e):
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    ce = np.maximum.accumulate(e)
    new = np.concatenate([[True], s[1:] > ce[:-1]])
    idx = np.flatnonzero(new)
    return s[idx], np.maximum.reduceat(e, idx)


def sweep(a_s, a_e, b_s, b_e):
    """intersection and difference (a minus b) of two sets of sorted, non-overlapping intervals"""
    pos = np.concatenate([a_s, a_e, b_s, b_e])
    da = np.concatenate([np.ones(len(a_s)), -np.ones(len(a_e)), np.zeros(len(b_s) + len(b_e))])
    db = np.concatenate([np.zeros(len(a_s) + len(a_e)), np.ones(len(b_s)), -np.ones(len(b_e))])
    o = np.argsort(pos, kind="stable")
    pos, ca, cb = pos[o], np.cumsum(da[o]), np.cumsum(db[o])
    s, e, ca, cb = pos[:-1], pos[1:], ca[:-1], cb[:-1]
    ok = e > s
    return (s[(ca > 0) & (cb > 0) & ok], e[(ca > 0) & (cb > 0) & ok]), (s[(ca > 0) & (cb == 0) & ok], e[(ca > 0) & (cb == 0) & ok])


# 1 kb windows with their GC class (windows with more than 100 N are dropped)
gcb = pd.read_csv(os.path.join(CL, "gc_bins_1kb.tsv"), sep="	")
gcb = gcb[(gcb["n"] <= 100) & (gcb["length"] == 1000)]
N_CLASS = 40
gcb["cls"] = np.minimum((gcb["gc"] / gcb["length"] * N_CLASS).astype(int), N_CLASS - 1)
gc_lookup = {}
for c, g in gcb.groupby("chrom"):
    arr = np.full(sizes[c] // 1000 + 2, -1, np.int16)
    arr[g["bin"].to_numpy()] = g["cls"].to_numpy()
    for s0, e0 in zip(*bl.get(c, EMPTY)):
        arr[s0 // 1000: e0 // 1000 + 1] = -1
    gc_lookup[cidx[c]] = arr
gc_pool = {}
for k, g in gcb.groupby("cls"):
    ch = np.array([cidx[c] for c in g["chrom"]], np.int64)
    bn = g["bin"].to_numpy(np.int64)
    keepm = np.ones(len(bn), bool)
    for c in np.unique(ch):
        sel = ch == c
        keepm[sel] = gc_lookup[c][bn[sel]] >= 0      # windows overlapping the blacklist were set to -1 above
    gc_pool[k] = (ch[keepm], bn[keepm])


def site_gc_class(keys):
    cls = np.full(len(keys), -1, np.int64)
    ci = keys // SHIFT
    for c, arr in gc_lookup.items():
        sel = ci == c
        cls[sel] = arr[(keys[sel] % SHIFT) // 1000]
    return cls


def place_gc(classes, width, rng):
    out = []
    for k, cnt in zip(*np.unique(classes[classes >= 0], return_counts=True)):
        chrs, bins = gc_pool[k]
        j = rng.integers(0, len(bins), size=cnt)
        out.append(chrs[j] * SHIFT + bins[j] * 1000 + rng.integers(0, 1000 - width, size=cnt))
    return np.concatenate(out)


# allowed placement segments: (chromosome code, start, end) for the whole genome, and inside / outside repeats
seg = {"all": [], "rep": [], "non": []}
for c in chroms:
    g_s, g_e = gaps.get(c, EMPTY)
    b_s, b_e = bl.get(c, EMPTY)
    excl = merge(np.concatenate([g_s, b_s]), np.concatenate([g_e, b_e])) if len(g_s) + len(b_s) else EMPTY
    allowed = sweep(np.array([0]), np.array([sizes[c]]), *excl)[1]
    r = reps.get(c, (np.array([], np.int64), np.array([], np.int64)))
    inside, outside = sweep(allowed[0], allowed[1], r[0], r[1])
    for k, v in (("all", allowed), ("rep", inside), ("non", outside)):
        seg[k].append((np.full(len(v[0]), cidx[c], np.int64), v[0], v[1]))
seg = {k: tuple(np.concatenate(x) for x in zip(*v)) for k, v in seg.items()}
seg_cum = {k: np.cumsum(v[2] - v[1]) for k, v in seg.items()}
print({k: int(seg_cum[k][-1]) for k in seg}, flush=True)


def place(kind, n, width, rng):
    ci, s, e = seg[kind]
    cum = seg_cum[kind]
    u = rng.integers(0, cum[-1], size=n)
    j = np.searchsorted(cum, u, side="right")
    off = u - (cum[j] - (e[j] - s[j]))
    return ci[j] * SHIFT + s[j] + off          # composite key chromosome * 2^32 + start


def downstream_distance(keys, width):
    keys = np.sort(keys)
    end = keys + width
    j = np.searchsorted(keys, end, side="right")
    ok = j < len(keys)
    jj = np.minimum(j, len(keys) - 1)
    ok &= (keys[jj] // SHIFT) == (keys // SHIFT)
    return (keys[jj] - end)[ok]


def in_set(keys, width, iv):
    out = np.zeros(len(keys), bool)
    for ci, c in enumerate(chroms):
        sel = (keys // SHIFT) == ci
        if not sel.any() or c not in iv:
            continue
        s, e = iv[c]
        pos = keys[sel] % SHIFT
        j = np.searchsorted(s, pos + width, side="left") - 1      # last repeat starting before the site end
        out[sel] = (j >= 0) & (e[np.maximum(j, 0)] > pos)
    return out


def cp(obs_log, ctrl_log):
    r = ks_2samp(obs_log, ctrl_log)
    return r.statistic * r.statistic_sign


inv = pd.read_csv(INV, encoding="utf-8-sig")
wanted = [m for m in args.motifs.split(",") if m] or list(inv["motif_id"])
tag = f"p{args.tier:.0e}".replace("-0", "-")
out_csv = os.path.join(CL, f"cluster_propensity_{tag}{'' if args.no_blacklist else '_bl'}_iter{args.iter}.csv")
rows = []
t0 = time.time()
for batch in sorted(glob.glob(os.path.join(FIMO, "*", "fimo.tsv")), key=lambda p: int(os.path.basename(os.path.dirname(p)).split("_")[0])):
    df = pd.read_csv(batch, sep="\t", comment="#", usecols=["motif_id", "sequence_name", "start", "stop", "p-value"],
                     dtype={"motif_id": str, "sequence_name": str, "start": np.int64, "stop": np.int64, "p-value": float})
    df = df[df["motif_id"].isin(wanted) & (df["p-value"] <= args.tier)]
    df["chrom"] = "chr" + df["sequence_name"]
    df = df[df["chrom"].isin(cidx)]
    for motif, g in df.groupby("motif_id", sort=False):
        rng = np.random.default_rng(abs(hash(motif)) % (2**32))
        width = int(round((g["stop"] - g["start"] + 1).median()))
        keys = np.array([cidx[c] for c in g["chrom"]], np.int64) * SHIFT + (g["start"].to_numpy() - 1)
        if bl:
            keys = keys[~in_set(keys, width, bl)]
        n = len(keys)
        d_obs = downstream_distance(keys, width)
        obs_log = np.log10(np.clip(d_obs, 1, None))
        rep_flag = in_set(keys, width, reps)
        n_rep = int(rep_flag.sum())
        gc_cls = site_gc_class(keys)
        cpa, cpb, cpc, medA, medB, medC = [], [], [], [], [], []
        for _ in range(args.iter):
            ka = place("all", n, width, rng)
            la = np.log10(np.clip(downstream_distance(ka, width), 1, None))
            cpa.append(cp(obs_log, la)); medA.append(np.median(la))
            kb = np.concatenate([place("rep", n_rep, width, rng), place("non", n - n_rep, width, rng)])
            lb = np.log10(np.clip(downstream_distance(kb, width), 1, None))
            cpb.append(cp(obs_log, lb)); medB.append(np.median(lb))
            kc = place_gc(gc_cls, width, rng)
            lc = np.log10(np.clip(downstream_distance(kc, width), 1, None))
            cpc.append(cp(obs_log, lc)); medC.append(np.median(lc))
        rows.append(dict(motif_id=motif, tier=args.tier, n_sites=n, width=width, frac_in_repeat=n_rep / n,
                         median_log10_dist_obs=np.median(obs_log), median_log10_dist_ctrl_paper=np.mean(medA),
                         median_log10_dist_ctrl_repeat=np.mean(medB), median_log10_dist_ctrl_gc=np.mean(medC),
                         mean_site_gc_class=float(gc_cls[gc_cls >= 0].mean() / N_CLASS), frac_within_20bp=float((d_obs <= 20).mean()),
                         frac_within_1kb=float((d_obs <= 1000).mean()),
                         CP_paper=np.mean(cpa), CP_paper_sd=np.std(cpa), CP_repeat_matched=np.mean(cpb), CP_repeat_matched_sd=np.std(cpb),
                         CP_gc_matched=np.mean(cpc), CP_gc_matched_sd=np.std(cpc)))
        pd.DataFrame(rows).to_csv(out_csv, index=False)
        print(f"{motif:24s} n={n:>9,}  CP paper={np.mean(cpa):+.3f}  CP repeat-matched={np.mean(cpb):+.3f}  CP GC-matched={np.mean(cpc):+.3f}  ({time.time() - t0:.0f}s)", flush=True)
print("done", out_csv)
