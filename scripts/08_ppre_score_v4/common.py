"""common.py -- shared paths, parameters, interval-overlap and statistics helpers for the PPARalpha LIVER / HEPATOCYTE ChIP-seq
validation of the PPRE (PPAR/RXR FIMO) score.  Imported by every numbered script; do not run directly.

The interval-overlap, DeLong-AUC, logistic-regression, decile and plotting helpers below are a copy of the tested helpers of
R_code/PPARA_validation_GSE244905/common.py (cardiomyocyte ChIP validation), kept in this folder so that this analysis is
self-contained and independent of a folder that is still being edited.

Conventions
-----------
* All genomic intervals are 0-based half-open [start0, end) (BED).  FIMO coordinates are 1-based inclusive (start0 = start - 1).
* Chromosomes are 'chr1'..'chr22','chrX','chrY' (UCSC style).  Internally a chromosome is an integer index 0..23 (X = 22, Y = 23;
  mouse uses 0..18, 22, 23) and an interval is encoded as one int64 (chrom_index << 32) + position so one global searchsorted works.
* ChIP-seq data are used strictly as an EXTERNAL validation set: nothing here changes, tunes or selects a PPRE score.
* Liver / hepatocyte ChIP-seq validates sequence-level PPAR/RXR binding potential only.  It says nothing about cardiomyocyte
  occupancy, regulatory activity, target-gene regulation or functional PPRE activity.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------------------------------- paths
SCRIPTS = Path(__file__).resolve().parent
PROJ = SCRIPTS.parent                                        # results/PPARA_liver_ChIP_validation
DATA = PROJ / "data"
TABLES = PROJ / "tables"
FIGS = PROJ / "figures"
REPORTS = PROJ / "reports"
CACHE = Path(os.environ.get("PPRE_CACHE_DIR", PROJ / "cache_local"))   # large intermediates; keep OFF OneDrive (set PPRE_CACHE_DIR)
for _d in (DATA, TABLES, FIGS, REPORTS, CACHE):
    _d.mkdir(parents=True, exist_ok=True)

ROOT = Path(r"C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship")
RCODE = ROOT / "R_code"
MOTIF_DIR = ROOT / "PPAR_RXR_motifs_Hocomoco_Jaspar"
FIMO_DIR = MOTIF_DIR / "FIMO_results_genome_masked"                       # raw FIMO (30 motifs, repeat-masked hg38, p <= 1e-4)
FIMO42_DIR = MOTIF_DIR / "FIMO_results_human42_unmasked"                  # 42 human-source motifs, unmasked hg38 (raw FIMO only)
OBJ1 = FIMO_DIR / "Objective1_analysis"
BED_DIR = OBJ1 / "Objective1_results" / "Bedfiles_output"
CHIPSEEKER_DIR = OBJ1 / "ChIPseeker_output"                                # per-motif annotated promoter hits (the scored hits)
HITS_BED = BED_DIR / "PPRE_FIMO_hits_30_motifs.bed"                        # genome-wide hits (15.9 M)
SCORE_BED = BED_DIR / "PPRE_CHIP_all_genes_test.bed"                       # final gene-level PPRE score table, 19,332 genes
MOTIF_FAMILY_CSV = MOTIF_DIR / "PPRE_threshold_sensitivity" / "results" / "motif_family_table.csv"
SCORING_XLSX = {
    "scoring_comparison": MOTIF_DIR / "Scoring_method_comparison" / "scoring_comparison_genome_masked.xlsx",
    "diversity_adjusted": MOTIF_DIR / "Scoring_method_comparison" / "diversity_adjusted_score" / "diversity_score_genome_masked.xlsx",
    "shrunk_average": MOTIF_DIR / "Scoring_method_comparison" / "shrunk_average" / "shrunk_average_genome_masked.xlsx",
}
RAW_PEAKS = DATA / "raw_peaks"

# --------------------------------------------------------------------------------------------- parameters
SEED = 20251005
WIN_UP, WIN_DOWN = 3000, 1000          # promoter = TSS-3000 .. TSS+1000 (strand-aware), as in the manuscript ChIPseeker call
DECAY = 3000                           # bp; decay constant of the manuscript score: score * exp(-|distanceToTSS| / 3000)
DECAY_BRIEF = 2000                     # bp; the value quoted in the analysis brief (sensitivity only)
STD_CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHR_IDX = {c: i for i, c in enumerate(STD_CHR)}

# --------------------------------------------------------------------------------------------- interval helpers
SHIFT = 32


def enc(chrom_idx, pos):
    """Encode (chromosome index, position) as one int64 so intervals of all chromosomes live on one number line."""
    return (np.asarray(chrom_idx, dtype=np.int64) << SHIFT) + np.asarray(pos, dtype=np.int64)


class IntervalSet:
    """Sorted, non-overlapping interval set (e.g. ChIP peaks) with vectorised overlap queries.

    `values` are optional per-interval numeric columns (e.g. peak signal) from which the maximum over all
    overlapping intervals is returned.
    """

    def __init__(self, chrom_idx, start0, end, **values):
        chrom_idx = np.asarray(chrom_idx); start0 = np.asarray(start0); end = np.asarray(end)
        o = np.lexsort((start0, chrom_idx))
        self.cidx, self.start, self.end = chrom_idx[o], start0[o], end[o]
        self.s = enc(self.cidx, self.start)
        self.e = enc(self.cidx, self.end)
        self.values = {k: np.asarray(v)[o] for k, v in values.items()}
        if len(self.s) > 1 and not (self.s[1:] >= self.e[:-1]).all():
            raise ValueError("intervals overlap each other; merge them first")

    def __len__(self):
        return len(self.s)

    def query(self, cidx, start0, end, with_values=()):
        """For each query interval: number of overlapping intervals (>= 1 bp) and the max of each requested value."""
        qs, qe = enc(cidx, start0), enc(cidx, end)
        lo = np.searchsorted(self.e, qs, side="right")      # first interval ending after the query start
        hi = np.searchsorted(self.s, qe, side="left")       # intervals starting before the query end
        n = np.maximum(hi - lo, 0)
        out = {"peak_idx": np.where(n > 0, lo, -1)}       # index of the first overlapping interval (cluster id of a hit)
        for name in with_values:
            v = self.values[name].astype(np.float64)
            m = np.full(len(qs), np.nan)
            for k in range(int(n.max()) if len(n) else 0):
                sel = n > k
                m[sel] = np.fmax(m[sel], v[lo[sel] + k])
            out[name] = m
        return n, out

    def overlaps(self, cidx, start0, end):
        return self.query(cidx, start0, end)[0] > 0

    def coverage(self, cidx, start0, end):
        """Number of query bases covered by this (non-overlapping) interval set."""
        if not hasattr(self, "_cum"):
            self._cum = np.concatenate([[0], np.cumsum(self.e - self.s)])

        def cov(x):
            k = np.searchsorted(self.s, x, side="left")          # intervals starting before x
            i = np.maximum(k - 1, 0)
            val = self._cum[i] + np.minimum(x, self.e[i]) - self.s[i]
            return np.where(k == 0, 0, val)
        return cov(enc(cidx, end)) - cov(enc(cidx, start0))


def merge_intervals(cidx, start0, end):
    """Merge overlapping (not merely touching) intervals; returns (cidx, start0, end) of the union."""
    cidx, start0, end = np.asarray(cidx), np.asarray(start0), np.asarray(end)
    o = np.lexsort((start0, cidx))
    cidx, start0, end = cidx[o], start0[o], end[o]
    s, e = enc(cidx, start0), enc(cidx, end)
    cm = np.maximum.accumulate(e)
    new = np.ones(len(s), bool)
    new[1:] = s[1:] >= cm[:-1]
    grp = np.cumsum(new) - 1
    first = np.flatnonzero(new)
    last_end = np.maximum.reduceat(end, first)
    return cidx[first], start0[first], last_end


# --------------------------------------------------------------------------------------------- statistics helpers
def mwu(pos, neg):
    """Mann-Whitney U (asymptotic, tie-corrected) with AUC = P(pos > neg) and Cliff's delta / rank-biserial = 2*AUC - 1."""
    r = stats.mannwhitneyu(pos, neg, alternative="two-sided", method="asymptotic")
    n1, n0 = len(pos), len(neg)
    auc = r.statistic / (n1 * n0)
    return dict(U=float(r.statistic), p=float(r.pvalue), auc=float(auc), cliffs_delta=float(2 * auc - 1))


def auc_delong(score, y, alpha=0.05):
    """ROC AUC with the exact DeLong variance (handles ties; fast for n0 >> n1 or n1 >> n0)."""
    score = np.asarray(score, dtype=np.float64); y = np.asarray(y).astype(bool)
    pos, neg = np.sort(score[y]), np.sort(score[~y])
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return dict(auc=np.nan, se=np.nan, lo=np.nan, hi=np.nan, n_pos=n1, n_neg=n0)
    # placement values: V10 for positives (share of negatives below, ties 0.5); V01 for negatives (share of positives above)
    v10 = (np.searchsorted(neg, pos, side="left") + np.searchsorted(neg, pos, side="right")) / (2.0 * n0)
    v01 = 1.0 - (np.searchsorted(pos, neg, side="left") + np.searchsorted(pos, neg, side="right")) / (2.0 * n1)
    auc = v10.mean()
    var = v10.var(ddof=1) / n1 + v01.var(ddof=1) / n0
    se = float(np.sqrt(var)); z = stats.norm.ppf(1 - alpha / 2)
    return dict(auc=float(auc), se=se, lo=float(auc - z * se), hi=float(auc + z * se), n_pos=n1, n_neg=n0)


def logit_fit(X, y, chunk=2_000_000, max_iter=60, tol=1e-9):
    """Logistic regression by Newton-Raphson (chunked, unpenalised).  X must contain an intercept column.
    Returns dict(beta, se, z, p, cov, loglik, converged)."""
    X = np.asarray(X); y = np.asarray(y, dtype=np.float64)
    n, k = X.shape
    beta = np.zeros(k)
    pbar = y.mean()
    beta[0] = np.log(pbar / (1 - pbar))
    converged = False
    for _ in range(max_iter):
        g = np.zeros(k); H = np.zeros((k, k))
        for i in range(0, n, chunk):
            xi = X[i:i + chunk].astype(np.float64); yi = y[i:i + chunk]
            p = 1.0 / (1.0 + np.exp(-(xi @ beta)))
            w = p * (1 - p)
            g += xi.T @ (yi - p)
            H += (xi * w[:, None]).T @ xi
        step = np.linalg.solve(H, g)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            converged = True
            break
    ll = 0.0
    H = np.zeros((k, k))
    for i in range(0, n, chunk):
        xi = X[i:i + chunk].astype(np.float64); yi = y[i:i + chunk]
        eta = xi @ beta
        ll += float(np.sum(yi * eta - np.logaddexp(0.0, eta)))
        p = 1.0 / (1.0 + np.exp(-eta)); w = p * (1 - p)
        H += (xi * w[:, None]).T @ xi
    cov = np.linalg.inv(H)
    se = np.sqrt(np.diag(cov)); z = beta / se
    return dict(beta=beta, se=se, z=z, p=2 * stats.norm.sf(np.abs(z)), cov=cov, loglik=ll, converged=converged)


def or_ci(beta, se, alpha=0.05):
    z = stats.norm.ppf(1 - alpha / 2)
    return float(np.exp(beta)), float(np.exp(beta - z * se)), float(np.exp(beta + z * se))


def logit_single(x, y, standardize=True):
    """Logistic regression of y on one (optionally standardised) predictor; returns coefficient, OR per SD, CI, p."""
    x = np.asarray(x, dtype=np.float64)
    mu, sd = x.mean(), x.std(ddof=1)
    z = (x - mu) / sd if standardize else x
    X = np.column_stack([np.ones(len(z)), z])
    f = logit_fit(X, y)
    orr, lo, hi = or_ci(f["beta"][1], f["se"][1])
    return dict(beta=float(f["beta"][1]), se=float(f["se"][1]), odds_ratio_per_SD=orr, or_ci_lo=lo, or_ci_hi=hi,
                p_value=float(f["p"][1]), sd=float(sd), mean=float(mu), converged=f["converged"])


def logit_fit_cluster(X, y, cluster, beta, chunk=2_000_000):
    """Cluster-robust (sandwich) covariance of a fitted logistic model.  `cluster` is an int array (length n);
    rows with cluster < 0 are their own cluster (independent), rows with the same id >= 0 share a cluster
    (here: all PPRE hits overlapping the same ChIP peak)."""
    X = np.asarray(X); y = np.asarray(y, dtype=np.float64); cluster = np.asarray(cluster)
    n, k = X.shape
    H = np.zeros((k, k)); Bs = np.zeros((k, k))
    gids = np.unique(cluster[cluster >= 0])
    gmap = {g: i for i, g in enumerate(gids)}
    G = np.zeros((len(gids), k))
    for i in range(0, n, chunk):
        xi = X[i:i + chunk].astype(np.float64); yi = y[i:i + chunk]; ci = cluster[i:i + chunk]
        p = 1.0 / (1.0 + np.exp(-(xi @ beta))); w = p * (1 - p)
        H += (xi * w[:, None]).T @ xi
        u = xi * (yi - p)[:, None]
        single = ci < 0
        Bs += u[single].T @ u[single]
        if (~single).any():
            idx = np.array([gmap[g] for g in ci[~single]])
            np.add.at(G, idx, u[~single])
    B = Bs + G.T @ G
    Hinv = np.linalg.inv(H)
    return Hinv @ B @ Hinv, len(gids)


def auc_cluster_boot(score, y, cluster, n_boot=2000, seed=SEED, alpha=0.05):
    """AUC with a percentile CI from a cluster bootstrap over the positives (clusters = ChIP peaks).  Negatives are treated as
    independent (their contribution to the variance is negligible when n_neg >> n_pos)."""
    score = np.asarray(score, dtype=np.float64); y = np.asarray(y).astype(bool); cluster = np.asarray(cluster)
    pos, neg = score[y], np.sort(score[~y])
    v10 = (np.searchsorted(neg, pos, side="left") + np.searchsorted(neg, pos, side="right")) / (2.0 * len(neg))
    cp = cluster[y]
    _, inv = np.unique(np.where(cp < 0, -1 - np.arange(len(cp)), cp), return_inverse=True)   # singletons get unique ids
    k = inv.max() + 1
    ssum = np.bincount(inv, weights=v10, minlength=k); msize = np.bincount(inv, minlength=k).astype(float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, k, size=(n_boot, k))
    boot = ssum[idx].sum(1) / msize[idx].sum(1)
    lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return dict(auc=float(v10.mean()), lo=float(lo), hi=float(hi), n_clusters=int(k))


def decile_boot(bins, y, cluster, nbins=10, n_boot=2000, seed=SEED):
    """Peak-cluster bootstrap CIs for the bin table: % positive per bin, odds ratio vs the lowest bin (0.5 continuity correction) and the
    Spearman correlation between bin number and % positive.  Positives are resampled by cluster (ChIP peak); bin sizes are fixed.
    Returns None when fewer than 5 clusters exist (a bootstrap is then meaningless)."""
    bins = np.asarray(bins); y = np.asarray(y).astype(bool); cluster = np.asarray(cluster)
    cp = cluster[y]
    ids = np.where(cp < 0, -1 - np.arange(len(cp)), cp)
    _, inv = np.unique(ids, return_inverse=True)
    G = int(inv.max()) + 1 if len(inv) else 0
    if G < 5:
        return None
    M = np.zeros((G, nbins)); np.add.at(M, (inv, bins[y] - 1), 1)
    n = np.bincount(bins, minlength=nbins + 1)[1:].astype(float)
    rng = np.random.default_rng(seed)
    C = np.empty((n_boot, nbins))
    for b0 in range(0, n_boot, 250):                       # chunked to bound memory
        idx = rng.integers(0, G, size=(min(250, n_boot - b0), G))
        C[b0:b0 + len(idx)] = M[idx].sum(1)
    pct = 100 * C / n
    odds = (C + 0.5) / (n - C + 0.5)
    orr = odds / odds[:, [0]]
    rk = stats.rankdata(pct, axis=1)
    x = np.arange(1, nbins + 1, dtype=float); x = (x - x.mean()) / x.std()
    rkc = (rk - rk.mean(1, keepdims=True)) / np.maximum(rk.std(1, keepdims=True), 1e-12)
    rho = (rkc * x).mean(1)
    q = lambda a: np.percentile(a, [2.5, 97.5], axis=0)
    return dict(pct=q(pct), orr=q(orr), rho=np.percentile(rho, [2.5, 97.5]))


def validate_score(score, y, cluster=None, nbins=10, seed=SEED, n_boot=2000):
    """Complete external validation of ONE continuous score against a binary ChIP-overlap label.

    Returns a dict with: descriptive stats per group, Mann-Whitney / Cliff's delta, AUC (DeLong CI and cluster-bootstrap CI),
    logistic regression on the standardised score (naive and cluster-robust SE), the decile table and its trend statistics.
    """
    score = np.asarray(score, dtype=np.float64); y = np.asarray(y).astype(bool)
    if cluster is None:
        cluster = -np.ones(len(y), dtype=np.int64)
    res = {}
    pos, neg = score[y], score[~y]
    res["n_pos"], res["n_neg"] = int(y.sum()), int((~y).sum())
    res["n_pos_clusters"] = int(len(np.unique(cluster[y & (cluster >= 0)])) + (y & (cluster < 0)).sum())
    res["desc_pos"], res["desc_neg"] = describe(pos), describe(neg)
    res["mwu"] = mwu(pos, neg)
    res["auc_delong"] = auc_delong(score, y)
    res["auc_cluster"] = auc_cluster_boot(score, y, cluster, n_boot=n_boot, seed=seed)
    mu, sd = score.mean(), score.std(ddof=1)
    X = np.column_stack([np.ones(len(score)), (score - mu) / sd])
    f = logit_fit(X, y)
    cov_r, ncl = logit_fit_cluster(X, y, cluster, f["beta"])
    se_r = float(np.sqrt(cov_r[1, 1])); z = stats.norm.ppf(0.975)
    b, se = float(f["beta"][1]), float(f["se"][1])
    res["logit"] = dict(beta=b, se_naive=se, odds_ratio_per_SD=float(np.exp(b)),
                        or_ci_lo_naive=float(np.exp(b - z * se)), or_ci_hi_naive=float(np.exp(b + z * se)),
                        p_naive=float(f["p"][1]), se_cluster=se_r,
                        or_ci_lo_cluster=float(np.exp(b - z * se_r)), or_ci_hi_cluster=float(np.exp(b + z * se_r)),
                        p_cluster=float(2 * stats.norm.sf(abs(b / se_r))), sd=float(sd), mean=float(mu), converged=f["converged"])
    bins = deciles_random_ties(score, nbins=nbins, seed=seed)
    tab = bin_table(bins, y, nbins=nbins)
    bt = decile_boot(bins, y, cluster, nbins=nbins, seed=seed)
    for col in ("pct_ci_lo_peakboot", "pct_ci_hi_peakboot", "or_ci_lo_peakboot", "or_ci_hi_peakboot"):
        tab[col] = np.nan
    res["trend"] = trend_stats(tab)
    res["trend"].update(spearman_ci_lo_peakboot=np.nan, spearman_ci_hi_peakboot=np.nan)
    if bt is not None:
        tab["pct_ci_lo_peakboot"], tab["pct_ci_hi_peakboot"] = bt["pct"][0], bt["pct"][1]
        tab["or_ci_lo_peakboot"], tab["or_ci_hi_peakboot"] = bt["orr"][0], bt["orr"][1]
        res["trend"].update(spearman_ci_lo_peakboot=float(bt["rho"][0]), spearman_ci_hi_peakboot=float(bt["rho"][1]))
    res["decile_table"] = tab
    res["bins"] = bins
    if res["n_pos_clusters"] < 5:                       # cluster-based intervals are not estimable with < 5 independent peaks
        res["auc_cluster"].update(lo=np.nan, hi=np.nan)
        for k_ in ("se_cluster", "or_ci_lo_cluster", "or_ci_hi_cluster", "p_cluster"):
            res["logit"][k_] = np.nan
    return res


def flatten_validation(res, label, extra=None):
    """One tidy row (dict) per validated score for the summary tables."""
    r = dict(analysis=label, **(extra or {}))
    r.update(n_total=res["n_pos"] + res["n_neg"], n_chip_overlap=res["n_pos"], n_independent_peaks=res["n_pos_clusters"],
             pct_chip_overlap=100 * res["n_pos"] / (res["n_pos"] + res["n_neg"]))
    for g in ("pos", "neg"):
        d = res["desc_" + g]
        tag = "ChIP_pos" if g == "pos" else "ChIP_neg"
        r.update({f"{tag}_median": d["median"], f"{tag}_mean": d["mean"], f"{tag}_q1": d["q1"], f"{tag}_q3": d["q3"], f"{tag}_IQR": d["iqr"]})
    r.update(mwu_U=res["mwu"]["U"], mwu_p=res["mwu"]["p"], cliffs_delta=res["mwu"]["cliffs_delta"])
    a, c = res["auc_delong"], res["auc_cluster"]
    r.update(auc=a["auc"], auc_ci_lo_delong=a["lo"], auc_ci_hi_delong=a["hi"], auc_ci_lo_peakboot=c["lo"], auc_ci_hi_peakboot=c["hi"])
    l = res["logit"]
    r.update(beta_per_SD=l["beta"], odds_ratio_per_SD=l["odds_ratio_per_SD"], or_ci_lo_naive=l["or_ci_lo_naive"], or_ci_hi_naive=l["or_ci_hi_naive"],
             p_naive=l["p_naive"], or_ci_lo_peakcluster=l["or_ci_lo_cluster"], or_ci_hi_peakcluster=l["or_ci_hi_cluster"],
             p_peakcluster=l["p_cluster"], score_SD=l["sd"])
    r.update(res["trend"])
    return r


def woolf_or(a, b, c, d, alpha=0.05):
    """Odds ratio (a/b)/(c/d) with Woolf CI; Haldane-Anscombe 0.5 correction when any cell is zero."""
    a, b, c, d = map(float, (a, b, c, d))
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    lor = np.log(a * d / (b * c)); se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = stats.norm.ppf(1 - alpha / 2)
    return float(np.exp(lor)), float(np.exp(lor - z * se)), float(np.exp(lor + z * se))


def wilson(k, n, alpha=0.05):
    if n == 0:
        return (np.nan, np.nan)
    z = stats.norm.ppf(1 - alpha / 2); p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(c - h), float(c + h)


def deciles_random_ties(score, nbins=10, seed=SEED):
    """Assign equal-sized bins (1 = lowest score) after ranking; tied scores are ordered by a seeded random key so that
    bin sizes are exactly equal and ties are not broken by file order (which would correlate with motif)."""
    n = len(score)
    rng = np.random.default_rng(seed)
    tiebreak = rng.random(n, dtype=np.float32)
    order = np.lexsort((tiebreak, score))
    rank = np.empty(n, dtype=np.int64); rank[order] = np.arange(n)
    return (rank * nbins // n + 1).astype(np.int8)


def bin_table(bins, y, nbins=10, label="decile"):
    """Per-bin counts, % positive (Wilson CI), enrichment and odds ratio (Woolf CI) relative to the lowest bin."""
    y = np.asarray(y).astype(bool)
    n = np.bincount(bins, minlength=nbins + 1)[1:]
    k = np.bincount(bins[y], minlength=nbins + 1)[1:]
    rows = []
    for i in range(nbins):
        lo, hi = wilson(k[i], n[i])
        rows.append({label: i + 1, "n": int(n[i]), "n_chip_overlap": int(k[i]), "pct_chip_overlap": 100 * k[i] / n[i],
                     "pct_ci_lo": 100 * lo, "pct_ci_hi": 100 * hi})
    t = pd.DataFrame(rows)
    base_pct = t.loc[0, "pct_chip_overlap"]
    t["enrichment_vs_lowest"] = t["pct_chip_overlap"] / base_pct if base_pct > 0 else np.nan
    ors = [woolf_or(t.loc[i, "n_chip_overlap"], t.loc[i, "n"] - t.loc[i, "n_chip_overlap"],
                    t.loc[0, "n_chip_overlap"], t.loc[0, "n"] - t.loc[0, "n_chip_overlap"]) for i in range(nbins)]
    t["odds_ratio_vs_lowest"] = [o[0] for o in ors]
    t["or_ci_lo"] = [o[1] for o in ors]
    t["or_ci_hi"] = [o[2] for o in ors]
    return t


def trend_stats(t, label="decile"):
    """Spearman correlation between bin number and % overlap (as requested) plus the Cochran-Armitage trend test."""
    rho, p = stats.spearmanr(t[label], t["pct_chip_overlap"])
    n = t["n"].to_numpy(float); k = t["n_chip_overlap"].to_numpy(float); x = t[label].to_numpy(float)
    N, K = n.sum(), k.sum(); pbar = K / N
    T = np.sum(x * (k - n * pbar))
    var = pbar * (1 - pbar) * (np.sum(n * x ** 2) - np.sum(n * x) ** 2 / N)
    z = T / np.sqrt(var) if var > 0 else np.nan
    return dict(spearman_rho=float(rho), spearman_p=float(p), cochran_armitage_z=float(z),
                cochran_armitage_p=float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan)


def describe(x):
    x = np.asarray(x, dtype=np.float64)
    q1, med, q3 = np.nanpercentile(x, [25, 50, 75])
    return dict(n=int(np.isfinite(x).sum()), median=float(med), mean=float(np.nanmean(x)), q1=float(q1), q3=float(q3),
                iqr=float(q3 - q1), sd=float(np.nanstd(x, ddof=1)))


# --------------------------------------------------------------------------------------------- plotting style
def set_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42,           # editable text in Illustrator
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7, "xtick.major.size": 3, "ytick.major.size": 3,
        "figure.dpi": 150, "savefig.dpi": 300, "axes.grid": False, "legend.frameon": False,
    })


# Okabe-Ito colour-blind-safe palette (consistent across all figures)
C_POS, C_NEG, C_MAIN, C_ALT, C_GREY = "#D55E00", "#56B4E9", "#0072B2", "#E69F00", "#7F7F7F"
C_GREEN, C_PURPLE, C_PINK = "#009E73", "#8B6BB1", "#CC79A7"


def save_fig(fig, name):
    """Save as vector PDF and 300-dpi PNG into figures/."""
    fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS / f"{name}.png", bbox_inches="tight", dpi=300)
    import matplotlib.pyplot as plt
    plt.close(fig)


def md_table(df, floatfmt="{:.4g}"):
    """Minimal GitHub-markdown table (avoids the optional `tabulate` dependency)."""
    def f(v):
        if isinstance(v, (float, np.floating)):
            return "NA" if not np.isfinite(v) else floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join(f(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return chr(10).join([head, sep] + body)


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return "<1e-300" if p < 1e-300 else (f"{p:.2e}" if p < 1e-3 else f"{p:.3f}")
