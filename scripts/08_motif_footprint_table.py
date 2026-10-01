"""Footprint table for the 42 human-source PPAR/RXR motifs.

For every motif (an "annotation unit" in the style of an annotation-footprint table):
  * total region length: union of all FIMO hit intervals genome-wide (overlapping and touching hits merged), in bp
  * footprint percentage: that length as a share of the scanned genome (sum of all sequence lengths in the FASTA)
  * unique genes: genes with at least one hit in the promoter window (3 kb upstream to 1 kb downstream of the TSS),
    from the ChIPseeker annotation of script 03b (sameStrand = FALSE)
and a second, pooled table: union of all 42 motifs, by database and by receptor.

Usage: python scripts/08_motif_footprint_table.py   (paths below are the local folders of the analysis)
"""
import csv, glob, os
import numpy as np
import pandas as pd

ROOT = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar"
FIMO_DIR = os.path.join(ROOT, "FIMO_results_human42_unmasked")
INV = os.path.join(ROOT, "Motifs_update_2026-09-30", "motif_inventory_human_only.csv")
CHIP_DIR = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline_full_local\ChIPseeker_output"
OUT = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline"
GENOME_BP = 3_099_750_718          # num-residues reported by FIMO for Homo_sapiens.GRCh38.dna.primary_assembly.fa


def union_length(starts, stops):
    """Total length of the union of 1-based inclusive intervals; touching intervals are merged."""
    if len(starts) == 0:
        return 0
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], stops[order]
    running_end = np.maximum.accumulate(e)
    new_group = np.empty(len(s), dtype=bool)
    new_group[0] = True
    new_group[1:] = s[1:] > running_end[:-1] + 1
    idx = np.flatnonzero(new_group)
    group_end = np.maximum.reduceat(e, idx)
    return int((group_end - s[idx] + 1).sum())


per_motif = {}                              # motif -> bp
intervals = {}                              # motif -> {chromosome: (starts, stops)}, kept to merge any group of motifs exactly
for tsv in sorted(glob.glob(os.path.join(FIMO_DIR, "*", "fimo.tsv"))):
    df = pd.read_csv(tsv, sep="	", comment="#", usecols=["motif_id", "sequence_name", "start", "stop"],
                     dtype={"motif_id": str, "sequence_name": str, "start": np.int32, "stop": np.int32})
    for motif, sub in df.groupby("motif_id", sort=False):
        total, per_chrom = 0, {}
        for chrom, ch in sub.groupby("sequence_name", sort=False):
            s, e = ch["start"].to_numpy(), ch["stop"].to_numpy()
            total += union_length(s, e)
            per_chrom[chrom] = (s, e)
        per_motif[motif] = total
        intervals[motif] = per_chrom
    print("done", os.path.basename(os.path.dirname(tsv)), flush=True)


def group_union_bp(motifs):
    chroms = {c for m in motifs for c in intervals[m]}
    total = 0
    for c in chroms:
        parts = [intervals[m][c] for m in motifs if c in intervals[m]]
        total += union_length(np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts]))
    return total


genes = {}
for f in glob.glob(os.path.join(CHIP_DIR, "*.txt")):
    motif = os.path.basename(f)[:-4]
    g = pd.read_csv(f, sep=" ", quotechar='"', usecols=lambda c: c == "ENSEMBL")["ENSEMBL"].dropna()
    genes[motif] = set(g.unique())

inv = list(csv.DictReader(open(INV, encoding="utf-8-sig")))
group = lambda r: "PPAR::RXR heterodimers" if r["category"] == "PPAR:RXR heterodimer" else r["name"].upper()
rank_g = ["PPAR::RXR heterodimers", "PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG"]
rank_db = ["JASPAR", "HOCOMOCO", "CIS-BP"]
inv.sort(key=lambda r: (rank_g.index(group(r)), rank_db.index(r["database"].split()[0]), r["motif_id"]))
ids = [r["motif_id"] for r in inv]

rows = [dict(unit="WGS", region_bp=GENOME_BP, footprint_pct=100.0, unique_genes="whole genome", n_motifs="")]
for r in inv:
    m = r["motif_id"]
    rows.append(dict(unit=m, region_bp=per_motif[m], footprint_pct=100 * per_motif[m] / GENOME_BP,
                     unique_genes=len(genes[m]), n_motifs=1))

# pooled table: all motifs, then by database, then by receptor / heterodimer
pooled = [("All 42 motifs", ids)]
for db in rank_db:
    pooled.append((f"{db} motifs", [r["motif_id"] for r in inv if r["database"].split()[0] == db]))
for g in rank_g:
    pooled.append((g if g != "PPAR::RXR heterodimers" else g, [r["motif_id"] for r in inv if group(r) == g]))
prow = [dict(unit="WGS", region_bp=GENOME_BP, footprint_pct=100.0, unique_genes="whole genome", n_motifs="")]
for name, members in pooled:
    bp = group_union_bp(members)
    prow.append(dict(unit=name, region_bp=bp, footprint_pct=100 * bp / GENOME_BP,
                     unique_genes=len(set().union(*[genes[m] for m in members])), n_motifs=len(members)))

os.makedirs(OUT, exist_ok=True)
for fname, data in (("motif_footprint_table.csv", rows), ("motif_footprint_table_pooled.csv", prow)):
    with open(os.path.join(OUT, fname), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
        w.writeheader()
        w.writerows(data)
print("motifs:", len(per_motif))
for r in prow:
    print(r)
