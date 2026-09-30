"""Summarise the FIMO run on the 42 human-source PPAR/RXR motifs into small tables.

Raw fimo.tsv files (several GB in total) are not stored in the repository; this script condenses them
into per-motif counts. Run inside WSL (or anywhere) with the FIMO output folder as argument:

    python3 06_fimo_human42_summary.py ~/PPRE_FIMO/FIMO_results_human42_unmasked results/fimo_human42
"""
import csv, glob, os, sys, collections

src, dst = sys.argv[1], sys.argv[2]
os.makedirs(dst, exist_ok=True)
tiers = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]

stats = collections.OrderedDict()
for tsv in sorted(glob.glob(os.path.join(src, "*", "fimo.tsv")), key=lambda p: int(os.path.basename(os.path.dirname(p)).split("_")[0])):
    batch = os.path.basename(os.path.dirname(tsv))
    with open(tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        ix = {c: i for i, c in enumerate(header)}
        for row in rd:
            if not row or row[0].startswith("#") or len(row) < len(header) - 1:
                continue
            mid = row[ix["motif_id"]]
            try:
                p = float(row[ix["p-value"]])
            except ValueError:
                continue
            s = stats.setdefault(mid, dict(batch=batch, n_hits=0, plus=0, minus=0, min_p=1.0, n_q_lt_0_05=0, max_score=-1e9,
                                            **{f"n_p_le_{t:g}": 0 for t in tiers}))
            s["n_hits"] += 1
            s["plus" if row[ix["strand"]] == "+" else "minus"] += 1
            s["min_p"] = min(s["min_p"], p)
            s["max_score"] = max(s["max_score"], float(row[ix["score"]]))
            q = row[ix["q-value"]] if "q-value" in ix and len(row) > ix["q-value"] else ""
            if q not in ("", "NA") and float(q) < 0.05:
                s["n_q_lt_0_05"] += 1
            for t in tiers:
                if p <= t:
                    s[f"n_p_le_{t:g}"] += 1
    print("done", batch, flush=True)

with open(os.path.join(dst, "fimo_hits_per_motif.csv"), "w", newline="") as fh:
    cols = ["motif_id"] + list(next(iter(stats.values())).keys())
    w = csv.writer(fh); w.writerow(cols)
    for mid, s in stats.items():
        w.writerow([mid] + [f"{v:.3g}" if isinstance(v, float) else v for v in s.values()])
print("motifs summarised:", len(stats), "| total hits:", sum(s["n_hits"] for s in stats.values()))
