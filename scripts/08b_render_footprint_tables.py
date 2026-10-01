"""Render the two footprint tables (per motif, and pooled) as figures in the style of an annotation-footprint table."""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\brouw\PPRE_cistrome_atlas\results\human42_pipeline"
plt.rcParams["font.family"] = ["Arial", "Liberation Sans", "DejaVu Sans"]


def nl_int(x):
    return f"{int(x):,}".replace(",", ".")


def nl_pct(x):
    x = float(x)
    s = f"{x:.2f}" if x >= 0.1 else f"{x:.3f}"
    return s.replace(".", ",")


def render(csv_name, out_name, pooled):
    rows = list(csv.DictReader(open(os.path.join(OUT, csv_name), encoding="utf-8-sig")))
    header = ["Annotation\nunit", "Total region length\n[bp]", "Footprint\npercentage [%]", "Unique genes\n[no.]"]
    if pooled:
        header.insert(1, "Motifs\n[no.]")
    body = []
    for r in rows:
        line = [r["unit"]]
        if pooled:
            line.append(r["n_motifs"] if r["n_motifs"] else "")
        line += [nl_int(r["region_bp"]), nl_pct(r["footprint_pct"]),
                 r["unique_genes"] if r["unique_genes"] == "whole genome" else nl_int(r["unique_genes"])]
        body.append(line)
    # column x positions (inches) and alignment
    widths = [2.6, 0.9, 2.2, 1.9, 1.6] if pooled else [2.6, 2.2, 1.9, 1.6]
    xs, x = [], 0.15
    for w in widths:
        xs.append(x)
        x += w
    total_w = x
    row_h, head_h = 0.27, 0.62
    height = head_h + row_h * len(body) + 0.45
    fig = plt.figure(figsize=(total_w, height))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, total_w); ax.set_ylim(0, height); ax.axis("off")
    top = height - 0.12
    ax.hlines(top, 0.1, total_w - 0.1, color="black", lw=1.4)
    for j, h in enumerate(header):
        cx = xs[j] + widths[j] / 2 if j else xs[0] + widths[0] / 2
        ax.text(cx, top - head_h / 2 - 0.02, h, ha="center", va="center", fontsize=10, fontweight="bold", linespacing=1.2)
    ax.hlines(top - head_h, 0.1, total_w - 0.1, color="black", lw=0.9)
    y = top - head_h - row_h * 0.62
    for i, line in enumerate(body):
        bold = line[0] in ("WGS", "All 42 motifs")
        for j, cell in enumerate(line):
            if j == 0:
                ax.text(xs[0] + 0.05, y, cell, ha="left", va="center", fontsize=9.5, fontweight="bold" if bold else "normal")
            else:
                ax.text(xs[j] + widths[j] - 0.35, y, cell, ha="right", va="center", fontsize=9.5, fontweight="bold" if bold else "normal")
        y -= row_h
        if line[0] == "WGS" or (pooled and line[0] == "All 42 motifs"):
            ax.hlines(y + row_h * 0.5, 0.1, total_w - 0.1, color="grey", lw=0.5)
    ax.hlines(0.17, 0.1, total_w - 0.1, color="black", lw=1.4)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{out_name}.{ext}"), dpi=300)
    plt.close(fig)
    print("saved", out_name)


render("motif_footprint_table.csv", "motif_footprint_table_per_motif", pooled=False)
render("motif_footprint_table_pooled.csv", "motif_footprint_table_pooled", pooled=True)
