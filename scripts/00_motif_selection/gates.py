import csv, re, math, os, json
import numpy as np

D = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar\Motifs_update_2026-09-30"
inv = list(csv.DictReader(open(os.path.join(D, "motif_inventory.csv"), encoding="utf-8-sig")))

# ---- read matrices from combined MEME
M = {}
cur = None
for line in open(os.path.join(D, "PPAR_RXR_all_motifs_combined.meme"), encoding="utf-8"):
    if line.startswith("MOTIF "):
        cur = line.split()[1]; M[cur] = []
    elif cur and re.match(r"^\s*[0-9.eE+-]+(\s+[0-9.eE+-]+){3}\s*$", line):
        M[cur].append([float(x) for x in line.split()])
M = {k: np.array(v) for k, v in M.items()}
assert len(M) == 71, len(M)

def ic(p):
    q = np.clip(p, 1e-12, 1)
    return float((2 + (p * np.log2(q)).sum(1)).sum())

def rc(p): return p[::-1, ::-1]

# ---- gate 1
g1 = {}
for x in inv:
    p = M[x["motif_id"]]
    zeros = float((p == 0).mean())
    icv = ic(p)
    n = None
    m = re.search(r"([0-9]+) sites", x["notes"])
    if x["database"].startswith("JASPAR"):
        pass
    if m: n = int(m.group(1))
    g1[x["motif_id"]] = dict(ic_bits=round(icv, 2), zero_frac=round(zeros, 3), nsites=n)

# ---- gate 2: pairwise similarity
ids = [x["motif_id"] for x in inv]
def best_sim(a, b):
    wa, wb = len(a), len(b)
    need = max(6, int(math.ceil(0.7 * max(wa, wb))))
    best, bestmad = -2, None
    for bb in (b, rc(b)):
        for off in range(-(wb - 1), wa):
            s0, e0 = max(0, off), min(wa, off + wb)
            L = e0 - s0
            if L < need: continue
            x = a[s0:e0].ravel(); y = bb[s0 - off:e0 - off].ravel()
            if x.std() == 0 or y.std() == 0: continue
            c = float(np.corrcoef(x, y)[0, 1])
            if c > best:
                best, bestmad = c, float(np.abs(a[s0:e0] - bb[s0 - off:e0 - off]).mean())
    return (best, bestmad) if best > -2 else (0.0, 1.0)

n = len(ids)
S = np.eye(n); MAD = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        c, mad = best_sim(M[ids[i]], M[ids[j]])
        S[i, j] = S[j, i] = c; MAD[i, j] = MAD[j, i] = mad
with open(os.path.join(D, "motif_pairwise_similarity.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow([""] + ids)
    for i in range(n): w.writerow([ids[i]] + [round(v, 3) for v in S[i]])

def cluster(th):
    # average linkage, similarity threshold
    cl = [[i] for i in range(n)]
    while True:
        best, pair = -1, None
        for a in range(len(cl)):
            for b in range(a + 1, len(cl)):
                s = np.mean([S[i, j] for i in cl[a] for j in cl[b]])
                if s > best: best, pair = s, (a, b)
        if best < th or pair is None: break
        a, b = pair; cl[a] = cl[a] + cl[b]; del cl[b]
    lab = {}
    for k, c in enumerate(sorted(cl, key=lambda c: min(c)), 1):
        for i in c: lab[ids[i]] = k
    return lab
cl98, cl95, cl90 = cluster(0.98), cluster(0.95), cluster(0.90)

# exact duplicate groups
dup = {}
gid = 0
for i in range(n):
    for j in range(i + 1, n):
        if S[i, j] >= 0.995 and MAD[i, j] <= 0.02 and len(M[ids[i]]) == len(M[ids[j]]):
            gi, gj = dup.get(ids[i]), dup.get(ids[j])
            if gi is None and gj is None:
                gid += 1; dup[ids[i]] = dup[ids[j]] = gid
            elif gi is None: dup[ids[i]] = gj
            elif gj is None: dup[ids[j]] = gi
            elif gi != gj:
                for k, v in list(dup.items()):
                    if v == gj: dup[k] = gi

# ---- gate 3: architecture (half-site geometry)
hs_plus = np.zeros((6, 4))  # RGGTCA
for k, s in enumerate(["AG", "G", "G", "T", "C", "A"]):
    for b in s: hs_plus[k, "ACGT".index(b)] = 1 / len(s)

def hs_scores(p):
    w = len(p); out = []
    for i in range(w - 5):
        f = float((p[i:i + 6] * hs_plus).sum(1).mean() * 1.0)
        r = float((p[i:i + 6] * rc(hs_plus)).sum(1).mean())
        out.append((f, r))
    return out

def architecture(p, th=0.55):
    sc = hs_scores(p)
    sites = []  # (pos, strand, score)
    for strand, idx in (("+", 0), ("-", 1)):
        cand = sorted([(s[idx], i) for i, s in enumerate(sc) if s[idx] >= th], reverse=True)
        for s, i in cand:
            if all(abs(i - j) >= 6 for _, j, st in [(0, q[0], q[1]) for q in sites]) or not sites:
                sites.append((i, strand, round(s, 2)))
    # global non-overlap: keep best by score
    sites.sort(key=lambda t: -t[2]); kept = []
    for i, st, s in sites:
        if all(abs(i - j) >= 6 for j, _, _ in kept): kept.append((i, st, s))
    kept.sort()
    if len(kept) == 0: return "no clear half-site", 0
    if len(kept) == 1: return "single half-site", 1
    if len(kept) >= 3: return "complex (>=3 half-sites)", len(kept)
    # take the two best-scoring
    top = sorted(sorted(kept, key=lambda t: -t[2])[:2])
    (i1, s1, _), (i2, s2, _) = top
    gap = i2 - (i1 + 6)
    if s1 == s2: return f"DR{gap}" + ("" if s1 == "+" else " (rev)"), len(kept)
    if s1 == "+" and s2 == "-": return f"IR{gap}", len(kept)
    return f"ER{gap}", len(kept)

# ---- species class
def species(x):
    db, src = x["database"], x["source_species"]
    if db.startswith("JASPAR"):
        return ("mouse" if "Mus" in src else "human"), "JASPAR species field"
    if db.startswith("HOCOMOCO"):
        return ("mouse" if src.startswith("Mouse") else "human"), "HOCOMOCO original_motif species"
    s = x["source"]
    if "Yin" in s: return "human", "inferred: Yin 2017 used human proteins"
    if "Gerstein" in s: return "human", "inferred: ENCODE human cell lines (H1-hESC, HepG2)"
    if "Kulakovskiy" in s: return "human", "inferred: HOCOMOCO v11 human motif (H11MO)"
    if "Isakova" in s: return "human", "inferred: SMiLE-seq of human TFs"
    if "Mathelier" in s: return "human", "inferred: JASPAR MA0066.2 is Homo sapiens"
    if "Vorontsov" in s: return "human", "inferred: GHT-SELEX with human genomic DNA"
    if "Jolma" in s: return "unknown (human or mouse clone)", "Jolma studies used both human and mouse clones; CIS-BP does not record which"
    if "Matys" in s: return "unknown (Transfac, mixed)", "Transfac entries mix species; not recorded by CIS-BP"
    return "unknown", "not recorded"

# ---- assemble
hoc_grade = lambda x: (re.search(r"quality (\w)", x["notes"]) or [None, None])[1]
rows = []
for x in inv:
    mid = x["motif_id"]; g = g1[mid]
    arch, nhs = architecture(M[mid])
    sp, basis = species(x)
    rows.append(dict(x, source_species_class=sp, species_basis=basis, hocomoco_quality=hoc_grade(x) or "",
                     ic_bits=g["ic_bits"], zero_fraction=g["zero_frac"], sparse_flag="SPARSE" if g["zero_frac"] >= 0.25 else "",
                     architecture=arch, n_half_sites=nhs, exact_duplicate_group=dup.get(mid, ""),
                     cluster_r095=cl95[mid], cluster_r090=cl90[mid]))
# propagate species through identical matrices (CIS-BP unknowns only)
byid = {r["motif_id"]: r for r in rows}
grp = {}
for r in rows:
    if r["exact_duplicate_group"] != "": grp.setdefault(r["exact_duplicate_group"], []).append(r)
for g, rs in grp.items():
    known = {r["source_species_class"] for r in rs if not r["source_species_class"].startswith("unknown") and not r["database"].startswith("CIS")}
    ids_known = [r["motif_id"] for r in rs if not r["database"].startswith("CIS") ]
    for r in rs:
        if r["source_species_class"].startswith("unknown"):
            if len(known) == 1:
                r["source_species_class"] = next(iter(known)) + " (via identical matrix)"; r["species_basis"] = "identical matrix to " + ", ".join(ids_known)
            elif len(known) > 1:
                r["source_species_class"] = "ambiguous (identical matrix listed for human and mouse)"; r["species_basis"] = "identical matrix to " + ", ".join(ids_known)
for r in rows: r["cluster_r098"] = cl98[r["motif_id"]]
# representative per r0.95 cluster
rank = lambda r: (r["sparse_flag"] == "SPARSE", {"A": 0, "B": 1, "C": 2, "D": 3, "": 1.5}[r["hocomoco_quality"]], -r["ic_bits"])
by = {}
for r in rows: by.setdefault(r["cluster_r095"], []).append(r)
for c, rs in by.items():
    rep = sorted(rs, key=rank)[0]
    for r in rs: r["cluster_r095_representative"] = "YES" if r is rep else ""
by2 = {}
for r in rows: by2.setdefault(r["cluster_r098"], []).append(r)
for c, rs in by2.items():
    rep = sorted(rs, key=rank)[0]
    for r in rs: r["cluster_r098_representative"] = "YES" if r is rep else ""
with open(os.path.join(D, "motif_inventory_gates.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

import collections
print("clusters r>=0.98:", len(set(cl98.values())));print("clusters r>=0.95:", len(set(cl95.values())), "| r>=0.90:", len(set(cl90.values())), "| exact-duplicate groups:", len(set(dup.values())), "motifs in them:", len(dup))
print("sparse:", [(r["motif_id"], r["zero_fraction"]) for r in rows if r["sparse_flag"]])
print("species:", collections.Counter((r["database"].split()[0], r["source_species_class"]) for r in rows))
print("architecture:", collections.Counter(r["architecture"] for r in rows))
print("zero_fraction distribution:", sorted(r["zero_fraction"] for r in rows)[::4])
