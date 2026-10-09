"""Rebuilds the motif prioritisation page with a chart-led layout (run from 'Claude outputs')."""
import pandas as pd, html, re, math
V4 = "../results/PPRE_score_v4/"
D = "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Minor internship/PPAR_RXR_motifs_Hocomoco_Jaspar/Motifs_update_2026-09-30/"
hl = pd.read_csv(V4 + "v4_standard_hitlevel.tsv", sep="\t"); gl = pd.read_csv(V4 + "v4_standard_gold.tsv", sep="\t"); ws = pd.read_csv(V4 + "v4_standard_windows.tsv", sep="\t")
inv = pd.read_csv(D + "motif_inventory_human_only.csv", encoding="utf-8-sig").set_index("motif_id"); sm = pd.read_csv(V4 + "motif_summary_42.csv").set_index("motif")
old = open("PPRE_motif_prioritisation.html", encoding="utf8").read()
css_base = old[:old.index("</style>")]                      # keeps <title>, font link and base tokens
NL = chr(10)
m = re.search(r"<section>" + NL + r"  <h2>Earlier tests: response and literature genes</h2>.*?</section>", old, re.S)
if m:
    earlier = m.group(0).replace("<section>" + NL + "  <h2>Earlier tests: response and literature genes</h2>", "").rsplit("</section>", 1)[0]
else:   # page already rebuilt: take the content of the collapsible block
    earlier = re.search(r"<summary>Earlier tests[^<]*</summary>" + NL + r"  (.*?)" + NL + r"</details>", old, re.S).group(1)

NAMES = [("reliable-25", "Reliable 25"), ("PPAR-centric (9)", "PPAR-centric 9"), ("no RXR half-sites (20)", "No RXR half-sites 20"), ("S0 remove RXR half-sites only (28)", "Half-sites only removed 28"),
         ("DR1 (22)", "DR1 22"), ("DR1-Q (20)", "DR1 no sparse 20"), ("DR1-NR (5)", "DR1 one per cluster 5"), ("SOFT weights", "Soft weights (all kept)")]
prim = gl[gl.gold == "human activation (primary)"]
W, LEFT, RIGHT, ROW, TOP = 380, 150, 16, 26, 26

def chart(title, xmin, xmax, ticks, fmt, ref, rows, xlabel, log=False):
    H = TOP + ROW * len(rows) + 34
    def X(v):
        if log: v, a, b = math.log(v), math.log(xmin), math.log(xmax)
        else: a, b = xmin, xmax
        return LEFT + (v - a) / (b - a) * (W - LEFT - RIGHT)
    g = f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}" style="width:100%;height:auto">'
    for t in ticks:
        g += f'<line x1="{X(t):.1f}" x2="{X(t):.1f}" y1="{TOP-8}" y2="{H-30}" stroke="var(--line)"/><text x="{X(t):.1f}" y="{H-16}" text-anchor="middle" font-size="10" fill="var(--muted)">{fmt(t)}</text>'
    g += f'<line x1="{X(ref):.1f}" x2="{X(ref):.1f}" y1="{TOP-8}" y2="{H-30}" stroke="var(--ink)" stroke-width="1.3"/>'
    for i, (lab, lo, mid, hi, tone) in enumerate(rows):
        y = TOP + ROW * i + 8
        g += f'<text x="{LEFT-8}" y="{y+4}" text-anchor="end" font-size="11" fill="var(--ink)">{html.escape(lab)}</text>'
        if lo is not None:
            g += f'<line x1="{X(lo):.1f}" x2="{X(hi):.1f}" y1="{y}" y2="{y}" stroke="var(--{tone})" stroke-width="2.4" stroke-linecap="round"/>'
        g += f'<circle cx="{X(mid):.1f}" cy="{y}" r="4.2" fill="var(--{tone})"/>' if mid is not None else f'<text x="{X(ref)+8:.1f}" y="{y+4}" font-size="10" fill="var(--muted)">n/a at hit level</text>'
    g += f'<text x="{(LEFT+W-RIGHT)/2:.1f}" y="{H-2}" text-anchor="middle" font-size="10" fill="var(--muted)">{html.escape(xlabel)}</text></svg>'
    return g

rows_hit = []
for key, lab in NAMES:
    d = hl[hl.set == key]
    if d.empty: rows_hit.append((lab, None, None, None, "muted")); continue
    rows_hit.append((lab, d.mh_or.min(), float(d.mh_or.median()), d.mh_or.max(), "good"))
rows_gene = []
for key, lab in NAMES:
    r = prim[prim.set == key].iloc[0]; rows_gene.append((lab, r.d_lo, r.d_auroc, r.d_hi, "warn"))
c1 = chart("Hit level: odds ratio of a hit lying in a ChIP peak, kept vs removed motifs", 0.8, 6, [1, 2, 3, 5], lambda v: f"{v:g}", 1, rows_hit, "odds ratio, kept vs removed hits", log=True)
c2 = chart("Gene level: AUROC difference vs all 42 motifs for verified PPAR targets", -0.12, 0.12, [-0.1, -0.05, 0, 0.05, 0.1], lambda v: f"{v:+.2f}".replace("+0.00", "0"), 0, rows_gene, "AUROC difference vs all 42 motifs")
base = prim[prim.set == "all 42"].iloc[0]

rel = sm[sm.reliable_25]; s2 = list(rel.index[rel.type.isin(["PPAR", "heterodimer"])]); s1 = list(rel.index[~((rel.type == "RXR") & (rel.length_bp <= 14))])
dr1 = inv.index[inv.architecture == "DR1"]; dr1q = [m for m in dr1 if str(inv.loc[m, "sparse_flag"]) != "SPARSE"]; nr = ["MA0855.1", "MA1148.1", "PPARA_M10653_3.10", "RXRA_M06774_3.10", "RXRG_M05829_3.10"]
def chips(l): return "".join(f'<code class="chip">{html.escape(m)}</code>' for m in sorted(l))
def sets_block(title, n, l, note): return f'<div class="setcard"><h4>{title} <span class="muted">({n})</span></h4><p class="note">{note}</p><div class="chips">{chips(l)}</div></div>'
cards = (sets_block("PPAR-centric", len(s2), s2, "PPAR and heterodimer models only.") + sets_block("DR1, one per cluster", len(nr), nr, "One motif per similarity cluster, DR1 only.")
         + sets_block("DR1 without sparse", len(dr1q), dr1q, "Direct-repeat models, sparse matrices removed.") + sets_block("No RXR half-sites", len(s1), s1, "Reliable motifs minus short RXR half-site models."))
rec = "".join(f"<tr><td>{g.split()[0]}</td><td class='n'>{int(gl[(gl.gold == g) & (gl.set == 'all 42')].n_pos.iloc[0])}</td><td class='n'>{gl[(gl.gold == g) & (gl.set == 'all 42')].auroc.iloc[0]:.3f}</td><td class='n'>{gl[(gl.gold == g) & (gl.set == 'PPAR-centric (9)')].auroc.iloc[0]:.3f}</td></tr>" for g in ["PPARA (all evidence)", "PPARD (all evidence)", "PPARG (all evidence)"])
win = "".join(f"<tr><td>{r.window}</td><td class='n'>{r.auroc:.3f}</td><td class='n'>{r.matched_auroc:.3f}</td></tr>" for r in ws.itertuples())

css = css_base.replace("<title>PPRE Motif Prioritisation</title>", "<title>PPRE Motif Prioritisation</title>") + '''
.two{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:16px 18px;min-width:0}
.card h3{font-family:var(--display);font-size:1.15rem;margin:0 0 6px}
.card.yes{border-top:4px solid var(--good)} .card.no{border-top:4px solid var(--warn)}
.big{font-family:var(--mono);font-size:1.6rem;font-weight:500;margin:2px 0 4px}
.muted{color:var(--muted);font-weight:400}
.chart{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:14px 14px 8px;min-width:0}
.chart h3{font-family:var(--display);font-size:1.1rem;margin:0 0 2px}
.chart p{font-size:.88rem;color:var(--muted);margin:0 0 6px}
.setcard{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:12px 14px;min-width:0}
.setcard h4{margin:0 0 2px;font-size:.98rem}
.chips{display:flex;flex-wrap:wrap;gap:6px}
.chip{background:var(--tint);border-radius:4px;padding:2px 6px;font-size:.78rem}
details{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:0 16px}
details>summary{cursor:pointer;padding:12px 0;font-family:var(--display);font-weight:600;font-size:1.1rem}
details[open]>summary{border-bottom:1px solid var(--line);margin-bottom:12px}
details .tbl{margin-bottom:12px}
</style>'''

body = f'''
<main>
<header>
  <p class="eyebrow">PPRE project · motifs</p>
  <h1>Can we prioritise the 42 motifs?</h1>
  <p class="sub">Seven ways to use fewer or weighted motifs, tested at hit level and against 136 verified human PPAR targets.</p>
</header>

<section class="two">
  <div class="card yes"><h3>Hits: more reliable</h3><div class="big">1.1 – 4.7×</div>
    <p>Hits of the kept motifs lie in ChIP peaks more often than hits of the removed motifs, at equal GC, on datasets not used for the choice. Also for sets chosen without any outcome data. Only removing the RXR half-sites alone does not help (0.9 – 1.7×).</p></div>
  <div class="card no"><h3>Genes: no measurable gain</h3><div class="big">AUROC {base.auroc:.2f}</div>
    <p>Against verified PPAR targets no set differs from all 42 motifs beyond noise, and all 42 are at chance for PPARD targets. So the score stays on all 42 motifs.</p></div>
</section>

<section>
  <h2>Hits versus genes, set by set</h2>
  <p>The same seven ways, two levels. Left of the black line means no better than the reference.</p>
  <div class="two">
    <div class="chart"><h3>Hit level</h3><p>Hit lies in a ChIP peak, kept vs removed motifs. Dot: median, bar: range over 6 ChIP sets.</p>{c1}</div>
    <div class="chart"><h3>Gene level</h3><p>Verified PPAR targets. Dot: estimate, bar: 95% interval.</p>{c2}</div>
  </div>
  <p class="note">Hit level: odds ratio stratified by flank GC, 6 held-out human ChIP sets; every bar stays above 1 except the half-sites-only set, whose range reaches 0.9. Gene level: 136 PPARgene activation targets (AUROC for all 42: {base.auroc:.3f}, matched negatives {base.matched_auroc:.3f}); every interval includes 0. With 136 targets a gain below about 0.04 cannot be resolved.</p>
</section>

<section>
  <h2>Why the gene score does not follow</h2>
  <ul>
    <li><strong>Known targets are often regulated from outside the promoter.</strong> The score reads only −3 kb to +1 kb of the TSS.</li>
    <li><strong>The signal sits close to the TSS.</strong> The narrowest window scores best and wider ones dilute it (table below).</li>
    <li><strong>The benchmark is small.</strong> Verified targets are few, and some receptors are not captured at all by promoter motifs (PPARD).</li>
  </ul>
</section>

<section>
  <h2>Which motifs are in each set</h2>
  <div class="two">{cards}</div>
  <p class="note">The full table of all 42 with reliability, receptor and hit counts is on the motif list page.</p>
</section>

<details>
  <summary>More tables: receptors and window</summary>
  <div class="tbl"><table>
    <thead><tr><th>Receptor (all evidence)</th><th>Targets</th><th>All 42 motifs</th><th>PPAR-centric (9)</th></tr></thead>
    <tbody>{rec}</tbody>
  </table></div>
  <div class="tbl"><table>
    <thead><tr><th>Window around the TSS</th><th>AUROC</th><th>Matched negatives</th></tr></thead>
    <tbody>{win}</tbody>
  </table></div>
</details>

<details>
  <summary>Earlier tests: hepatocyte response, literature genes, cardiomyocytes</summary>
  {earlier}
</details>

<section>
  <h2>What to say in the paper</h2>
  <ul>
    <li>The score uses all 42 motifs, and all hits are released. Each hit can carry a flag for the reliable subset, since reliable-subset hits are clearly more often bound.</li>
    <li>Reliable and DR1 subsets are a sensitivity analysis: they raise binding on held-out data without a measurable gain for gene-level prediction.</li>
    <li>Soft weights are an option that keeps every motif and every hit.</li>
    <li>The strongest claim is site prioritisation at hit level. Gene-level scoring is the weak part.</li>
  </ul>
</section>

<footer>Details and pre-set rules: <code>results/PPRE_score_v4/DESIGN_V4.md</code> (Addenda 5, 10, 11 to 15).</footer>
</main>
'''
open("PPRE_motif_prioritisation.html", "w", encoding="utf8").write(css + body)
print("built", len(css + body))
