import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

D = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar\Motifs_update_2026-09-30"
inv = {r["motif_id"]: r for r in csv.DictReader(open(D + r"\motif_inventory_human_only.csv", encoding="utf-8-sig"))}
rows = list(csv.reader(open(D + r"\motif_pairwise_similarity.csv")))
allids = rows[0][1:]
S = np.array([[float(v) for v in r[1:]] for r in rows[1:]])
ids = [i for i in allids if i in inv]
ix = [allids.index(i) for i in ids]
S = S[np.ix_(ix, ix)]
np.fill_diagonal(S, 1)
Z = linkage(squareform(1 - S, checks=False), method="average")
order = dendrogram(Z, no_plot=True)["leaves"]
c98 = fcluster(Z, 0.02, "distance"); c95 = fcluster(Z, 0.05, "distance")
print("clusters r>=0.98:", len(set(c98)), "| r>=0.95:", len(set(c95)))
oid = [ids[i] for i in order]
So = S[np.ix_(order, order)]

dbcol = {"JASPAR": "#4C78A8", "HOCOMOCO": "#F58518", "CIS-BP": "#54A24B"}
archcol = {"DR1": "#B279A2", "IR0": "#E45756", "IR3": "#9D755D", "single half-site": "#BAB0AC", "complex (>=3 half-sites)": "#FFD92F"}
def db(i): return inv[i]["database"].split()[0]

fig = plt.figure(figsize=(13, 11))
gs = fig.add_gridspec(2, 4, width_ratios=[1.3, 0.12, 0.12, 6], height_ratios=[0.1, 1], wspace=0.03, hspace=0.02)
axd = fig.add_subplot(gs[1, 0]); axa = fig.add_subplot(gs[1, 1]); axb = fig.add_subplot(gs[1, 2]); axh = fig.add_subplot(gs[1, 3])
dendrogram(Z, orientation="left", ax=axd, no_labels=True, color_threshold=0, above_threshold_color="#555", link_color_func=lambda k: "#555")
axd.invert_yaxis()
for th, ls in ((0.02, "--"), (0.05, ":")):
    axd.axvline(th, color="crimson", ls=ls, lw=1)
axd.set_xlim(0.45, 0); axd.set_xlabel("1 − r"); axd.set_yticks([])
for s in ("top", "right", "left"): axd.spines[s].set_visible(False)
axd.set_ylim(len(oid) * 10, 0)
for ax, cols, name in ((axa, [archcol[inv[i]["architecture"]] for i in oid], "arch"), (axb, [dbcol[db(i)] for i in oid], "db")):
    for k, c in enumerate(cols): ax.add_patch(plt.Rectangle((0, k), 1, 1, color=c))
    ax.set_xlim(0, 1); ax.set_ylim(len(oid), 0); ax.axis("off")
im = axh.imshow(So, cmap="viridis", vmin=0.7, vmax=1.0, aspect="auto", interpolation="nearest")
axh.set_xticks(range(len(oid))); axh.set_yticks(range(len(oid)))
lab = [f"{i}{'  [sparse]' if inv[i]['sparse_flag'] else ''}" for i in oid]
axh.set_xticklabels(lab, rotation=90, fontsize=7); axh.yaxis.tick_right(); axh.set_yticklabels(lab, fontsize=7)
# outline r>=0.98 clusters (contiguous in dendrogram order)
lab98 = [c98[order[k]] for k in range(len(oid))]
start = 0
for k in range(1, len(oid) + 1):
    if k == len(oid) or lab98[k] != lab98[start]:
        if k - start > 1:
            axh.add_patch(plt.Rectangle((start - .5, start - .5), k - start, k - start, fill=False, ec="white", lw=1.6))
        start = k
cb = fig.colorbar(im, ax=axh, fraction=0.025, pad=0.16); cb.set_label("best-alignment Pearson r (values < 0.7 shown as 0.7)")
handles = [Patch(color=v, label=k) for k, v in archcol.items()] + [Patch(color=v, label=k) for k, v in dbcol.items()]
fig.legend(handles=handles, loc="upper left", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.02, 0.97), title="Left bars: architecture, then database", title_fontsize=8)
fig.suptitle(f"Gate 2: motif redundancy, 42 human-source PPAR/RXR motifs\nAverage-linkage clustering; dashed/dotted = r 0.98 / 0.95 cut ({len(set(c98))} / {len(set(c95))} clusters); white boxes = r≥0.98 clusters", y=0.995, fontsize=11)
fig.subplots_adjust(top=0.90, bottom=0.2, left=0.02, right=0.93)
for ext in ("png", "pdf"): fig.savefig(D + rf"\gate2_redundancy_heatmap_human42.{ext}", dpi=170)
print("saved")
