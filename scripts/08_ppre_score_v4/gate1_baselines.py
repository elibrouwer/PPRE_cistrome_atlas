"""v4 Gate 1: baselines against human hepatocyte agonist responses (GSE17251 Wy14643, GSE53399 GW7647).
Labels fixed before scoring: up = top 5% of paired t (agonist vs vehicle, same donor). Exploratory development sets, not the locked set."""
import gzip, re, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
rng = np.random.default_rng(20261008)
R = "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/R_code"
G0 = "C:/Users/brouw/AppData/Local/Temp/claude/C--Users-brouw-OneDrive---Universiteit-Utrecht-Master-Bioinformatics-Minor-internship-R-code/f7175609-6765-498e-8602-f582703f131b/scratchpad/gate0/"
C = "C:/Users/brouw/AppData/Local/Temp/claude/C--Users-brouw-OneDrive---Universiteit-Utrecht-Master-Bioinformatics-Minor-internship-R-code/83a53ee4-67d0-4044-980d-5b615b596385/scratchpad/cache/human/"
OUT = f"{R}/results/PPRE_score_v4/"

def matrix(path):
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    titles = [t.strip('"') for t in next(l for l in lines if l.startswith("!Sample_title")).split("\t")[1:]]
    i = next(k for k, l in enumerate(lines) if l.startswith('"ID_REF"'))
    rows = [l.split("\t") for l in lines[i + 1:] if not l.startswith("!")]
    m = pd.DataFrame([r[1:] for r in rows], index=[r[0].strip('"') for r in rows], columns=titles).astype(float)
    return m

_n = next(i for i, l in enumerate(gzip.open(G0 + "GPL13158.annot.gz", "rt", errors="replace")) if l.startswith("ID\t"))
ann = pd.read_csv(G0 + "GPL13158.annot.gz", sep="\t", skiprows=_n, low_memory=False, usecols=["ID", "Gene symbol"], quoting=3, on_bad_lines="skip")
ann["sym"] = ann["Gene symbol"].astype(str).str.split("///").str[0]
ann = ann[(ann.sym != "nan") & (ann.sym != "")]
ann["ID"] = ann.ID.str.replace("_PM", "", regex=False)
pm = ann.drop_duplicates("ID").set_index("ID").sym
def to_gene(m):
    m = m.copy(); m.index = [i.replace("_PM", "") for i in m.index]
    m = m.loc[m.index.intersection(pm.index)].copy()
    m["sym"] = pm.loc[m.index].values
    m["mu"] = m.drop(columns="sym").mean(axis=1)
    m = m.sort_values("mu", ascending=False).drop_duplicates("sym")
    return m.set_index("sym").drop(columns="mu")

def paired_t(m, pairs):
    a = m[[p[0] for p in pairs]].to_numpy(); b = m[[p[1] for p in pairs]].to_numpy()
    d = a - b; t, p = stats.ttest_1samp(d, 0, axis=1)
    return pd.DataFrame({"stat": t, "lfc": d.mean(1), "expr": m.mean(axis=1).to_numpy()}, index=m.index)

# GSE17251: 24h Wy14643 vs DMSO, 6 donors
m1 = to_gene(matrix(G0 + "GSE17251_series_matrix.txt.gz"))
pairs1 = [(f"Donor{d}_24h_Wy", f"Donor{d}_24h_DMSO") for d in range(1, 7)]
lab1 = paired_t(m1, pairs1)
# GSE53399: GW7647 10uM vs 0uM at 24h, donors with both (names Hu####_24hr_dose, no T-prefix)
m2raw = matrix(G0 + "GSE53399_series_matrix.txt.gz")
cols = [c for c in m2raw.columns if re.fullmatch(r"Hu\d+_24hr_[\d.]+uM", c)]
donors = sorted({c.split("_")[0] for c in cols if f"{c.split('_')[0]}_24hr_10uM" in cols and f"{c.split('_')[0]}_24hr_0uM" in cols})
print("GSE53399 donors with 24h 10uM and 0uM:", donors)
m2 = to_gene(m2raw[cols])
lab2 = paired_t(m2, [(f"{d}_24hr_10uM", f"{d}_24hr_0uM") for d in donors])
lab1.to_csv(OUT + "labels_GSE17251.csv"); lab2.to_csv(OUT + "labels_GSE53399.csv")
sets = {"GSE17251 hepatocyte Wy14643 24h": lab1, "GSE53399 hepatocyte GW7647 24h 10uM": lab2}

# scores
sc = pd.read_csv(f"{R}/results/PPRE_score_v3/tables/v3_gene_scores_all_variants_42.csv")
gi = pd.read_parquet(C + "genes.parquet")[["ENSEMBL", "hit_count"]]
tss = pd.read_parquet(C + "tss.parquet")
chrom = {i: f"chr{c}" for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
# promoter ATAC: heart LV ATAC signal (ENCFF174HPZ), promoter -3000..+1000 (as v3)
atac = pd.read_csv(f"{R}/benchmark/external_data/heart_chromatin/ENCFF174HPZ.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 6])
atac.columns = ["c", "s", "e", "sig"]
res = []
by = {c: g.sort_values("s") for c, g in atac.groupby("c")}
tss["ch"] = tss.cidx.map(chrom)
tss["ps"] = np.where(tss.strand == "+", tss.tss - 3000, tss.tss - 1000); tss["pe"] = np.where(tss.strand == "+", tss.tss + 1000, tss.tss + 3000)
vals = np.zeros(len(tss))
for c, g in tss.groupby("ch"):
    a = by.get(c)
    if a is None: continue
    s, e, sig = a.s.to_numpy(), a.e.to_numpy(), a.sig.to_numpy()
    for i, (ps, pe) in zip(g.index, zip(g.ps, g.pe)):
        k = np.searchsorted(s, pe)
        j0 = max(0, np.searchsorted(s, ps - 5000))  # peaks are short
        ss, ee, gg = s[j0:k], e[j0:k], sig[j0:k]
        ov = np.clip(np.minimum(ee, pe) - np.maximum(ss, ps), 0, None)
        vals[tss.index.get_loc(i)] = (ov * gg).sum() / 4000
tss["atac"] = vals
atac_g = tss.groupby("ENSEMBL").atac.max().rename("atac_prom")
sc = sc.merge(gi, on="ENSEMBL", how="left").merge(atac_g, on="ENSEMBL", how="left")
sc["atac_prom"] = sc.atac_prom.fillna(0)
sc["log_atac"] = np.log1p(sc.atac_prom)
z = lambda a: (np.asarray(a, float) - np.nanmean(a)) / np.nanstd(a, ddof=1)
# simplest v4 candidate: v3 motif evidence plus promoter accessibility, equal-weight z-sum (no fitted weights)
sc["v3_plus_atac"] = z(sc.v3) + z(sc.log_atac)
sc["gc_only"] = sc.gc; sc["length_only"] = np.log10(sc.prom_len); sc["tss_count_only"] = sc.n_tss
sc["hit_count"] = sc.hit_count.fillna(0)
SC = {"v3": "v3", "v3-noTSS": "v3_noTSS", "original sum": "original_sum_30", "hit count": "hit_count", "GC only": "gc_only",
      "promoter length only": "length_only", "#TSS only": "tss_count_only", "ATAC only": "log_atac", "v3 + ATAC (v4 simplest)": "v3_plus_atac"}
u = sc[sc.SYMBOL.map(sc.SYMBOL.value_counts()) == 1].set_index("SYMBOL")

def auc_of(x, y):
    r = stats.rankdata(x); n1 = y.sum(); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (~y).sum())
rows = []
for sn, lab in sets.items():
    idx = lab.index.intersection(u.index); d = u.loc[idx].copy(); L = lab.loc[idx]
    n = len(d); import os; FR = float(os.environ.get("FRAC", ".05")); k = int(round(FR * n)); y = (L.stat.rank(ascending=False) <= k).to_numpy()
    cov = np.column_stack([z(np.log10(d.prom_len)), z(d.gc), z(np.log(d.n_tss)), z(L.expr)])
    pi, ni = np.flatnonzero(y), np.flatnonzero(~y); print(sn, "genes", n, "up", y.sum())
    for name, col in SC.items():
        x = d[col].to_numpy(float); x = np.where(np.isfinite(x), x, np.nanmedian(x))
        a = auc_of(x, y); bs = []
        for _ in range(1000):
            aa = rng.choice(pi, len(pi)); bb = rng.choice(ni, len(ni)); bs.append(auc_of(np.r_[x[aa], x[bb]], np.r_[np.ones(len(aa), bool), np.zeros(len(bb), bool)]))
        if col in ("gc_only", "length_only", "tss_count_only"):  # score is itself a covariate: adjusted OR undefined
            rows.append(dict(dataset=sn, score=name, auc=a, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), adj_or=np.nan, or_lo=np.nan, or_hi=np.nan)); continue
        f = sm.Logit(y.astype(int), sm.add_constant(np.column_stack([z(x), cov]))).fit(disp=0)
        rows.append(dict(dataset=sn, score=name, auc=a, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                         adj_or=np.exp(f.params[1]), or_lo=np.exp(f.params[1] - 1.96 * f.bse[1]), or_hi=np.exp(f.params[1] + 1.96 * f.bse[1])))
o = pd.DataFrame(rows); o.to_csv(OUT + f"gate1_baselines_frac{FR}.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(o.round(3).to_string(index=False))
sc.to_csv(OUT + "gate1_scores_with_atac.csv", index=False)
