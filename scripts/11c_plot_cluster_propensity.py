"""Dumbbell plot of the cluster propensity (CP) per motif under the three controls (paper, repeat-matched, GC-matched)."""
import glob, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CL = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline\clustering"
INV = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar\Motifs_update_2026-09-30\motif_inventory_human_only.csv"
inv = pd.read_csv(INV, encoding="utf-8-sig").set_index("motif_id")
files = {t: sorted(glob.glob(os.path.join(CL, f"cluster_propensity_{t}_bl_iter*.csv")) or glob.glob(os.path.join(CL, f"cluster_propensity_{t}_iter*.csv")))[-1] for t in ("p1e-4", "p1e-6")
         if glob.glob(os.path.join(CL, f"cluster_propensity_{t}_*iter*.csv"))}
files = {t: f for t, f in files.items() if len(pd.read_csv(f)) >= 40}
print("tiers plotted:", list(files))

order = pd.read_csv(os.path.join(CL, files[list(files)[0]])) if False else None
base = pd.read_csv(files[list(files)[0]]).copy()
base["family"] = ["heterodimer" if inv.loc[m, "category"] == "PPAR:RXR heterodimer" else ("PPAR half-site" if str(inv.loc[m, "name"]).upper().startswith("PPAR") else "RXR half-site") for m in base.motif_id]
base = base.sort_values("CP_paper")
motifs = list(base.motif_id)

fig, axes = plt.subplots(1, len(files), figsize=(5.6 * len(files) + 1.6, 11), sharey=True, squeeze=False)
col = {"paper": "#8a8a8a", "repeat": "#3B7FB6", "gc": "#C2561F"}
for ax, (tier, f) in zip(axes[0], files.items()):
    d = pd.read_csv(f).set_index("motif_id").reindex(motifs)
    y = range(len(d))
    ax.hlines(y, d.CP_gc_matched, d.CP_paper, color="#c9c2b4", lw=2, zorder=1)
    ax.scatter(d.CP_paper, y, s=34, color=col["paper"], label="paper control (random, outside gaps)", zorder=3)
    ax.scatter(d.CP_repeat_matched, y, s=26, facecolors="none", edgecolors=col["repeat"], label="repeat-matched", zorder=3)
    ax.scatter(d.CP_gc_matched, y, s=34, color=col["gc"], label="GC-matched", zorder=4)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title(f"p ≤ {tier[1:].replace('e-', 'e−')}  ({int(d.n_sites.sum()):,} sites)", fontsize=11, loc="left")
    ax.set_xlabel("cluster propensity (signed KS statistic; > 0 = sites closer than control)")
    ax.grid(axis="x", color="#e8e4da", lw=0.6)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0][0].set_yticks(list(range(len(motifs))))
axes[0][0].set_yticklabels([f"{m}  [{base.set_index('motif_id').loc[m, 'family'].split()[0]}]" for m in motifs], fontsize=8)
axes[0][0].legend(loc="lower right", fontsize=8, frameon=False)
fig.suptitle("Genomic clustering of PPAR/RXR motif sites, 42 human-source motifs", x=0.01, ha="left", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(CL, f"cluster_propensity_dumbbell.{ext}"), dpi=200)
print("saved")
