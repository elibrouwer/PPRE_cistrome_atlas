"""Adds the Addendum 15 benchmark results to the motif prioritisation page and the summary pointer (run from 'Claude outputs')."""
import pandas as pd, html
V4 = "../results/PPRE_score_v4/"
hl = pd.read_csv(V4 + "v4_standard_hitlevel.tsv", sep="\t"); gl = pd.read_csv(V4 + "v4_standard_gold.tsv", sep="\t"); ws = pd.read_csv(V4 + "v4_standard_windows.tsv", sep="\t")
order = ["reliable-25", "PPAR-centric (9)", "no RXR half-sites (20)", "S0 remove RXR half-sites only (28)", "DR1 (22)", "DR1-Q (20)", "DR1-NR (5)"]
lab = {"S0 remove RXR half-sites only (28)": "Only RXR half-sites removed (28)"}
rows = ""
for s in order:
    d = hl[hl.set == s]
    rows += f"<tr><td>{html.escape(lab.get(s, s))}</td><td class='n'>{d.mh_or.min():.1f}–{d.mh_or.max():.1f}</td><td class='n'>{int((d.lo > 1).sum())} of 6</td></tr>\n"
prim = gl[gl.gold == "human activation (primary)"]; base = prim[prim.set == "all 42"].iloc[0]
grow = ""
for s in ["PPAR-centric (9)", "DR1 (22)", "SOFT weights", "no RXR half-sites (20)", "reliable-25", "DR1-Q (20)", "S0 remove RXR half-sites only (28)", "DR1-NR (5)"]:
    r = prim[prim.set == s].iloc[0]; cls = "warn" if r.d_lo < 0 else "good"
    grow += f"<tr><td>{html.escape(lab.get(s, s))}</td><td class='n'>{r.auroc:.3f}</td><td class='n'>{r.d_auroc:+.3f} ({r.d_lo:+.3f} to {r.d_hi:+.3f})</td><td class='n'>{r.d_matched:+.3f}</td><td><span class='pill {cls}'>{'interval includes 0' if r.d_lo < 0 else 'above 0'}</span></td></tr>\n"
rec = ""
for g in ["PPARA (all evidence)", "PPARD (all evidence)", "PPARG (all evidence)"]:
    a = gl[(gl.gold == g) & (gl.set == "all 42")].iloc[0]; b = gl[(gl.gold == g) & (gl.set == "PPAR-centric (9)")].iloc[0]
    rec += f"<tr><td>{g.split()[0]}</td><td class='n'>{int(a.n_pos)}</td><td class='n'>{a.auroc:.3f}</td><td class='n'>{b.auroc:.3f}</td></tr>\n"
w = "".join(f"<tr><td>{r.window}</td><td class='n'>{r.auroc:.3f}</td><td class='n'>{r.matched_auroc:.3f}</td></tr>\n" for r in ws.itertuples())
sec = f'''<section>
  <h2>Tested the way other research does it</h2>
  <p>The first tests looked at gene response in two small datasets. Published motif and target benchmarks work differently, so we repeated the comparison with their practice: judge <strong>hits</strong> against bound and unbound sequences at equal GC, use a larger set of <strong>verified targets</strong> with precision-recall and matched negatives, and sweep the linking window. The rule was written before the run.</p>
  <h3>1. Are the kept hits more reliable? Yes.</h3>
  <p>For each ChIP dataset not used to choose the motifs, we compared hits of the kept motifs with hits of the removed ones at the same GC. A ratio above 1 means kept hits fall in ChIP peaks more often.</p>
  <div class="tbl"><table>
    <thead><tr><th>Motif set</th><th>Odds ratio kept vs removed (range over 6 ChIP sets)</th><th>Lower bound above 1</th></tr></thead>
    <tbody>
{rows}    </tbody>
  </table></div>
  <p class="note">Example: with the reliable 25, precision is 0.021 against 0.005 for removed hits in PPARG liver, and 0.102 against 0.069 in RXRA liver. This also holds for sets chosen without any outcome data (DR1, similarity clusters).</p>
  <h3>2. Does it improve gene-level prediction of verified targets? Not measurably.</h3>
  <p>Reference: 136 verified human PPAR targets with activation evidence from PPARgene. All 42 motifs give AUROC {base.auroc:.3f} (matched negatives {base.matched_auroc:.3f}) and a precision-recall score of {base.auprc:.3f} against a baseline of {base.auprc_baseline:.3f}.</p>
  <div class="tbl"><table>
    <thead><tr><th>Motif set</th><th>AUROC</th><th>Difference vs all 42 (95% interval)</th><th>Matched negatives</th><th>Verdict</th></tr></thead>
    <tbody>
{grow}    </tbody>
  </table></div>
  <div class="tbl" style="margin-top:12px"><table>
    <thead><tr><th>Receptor (all evidence)</th><th>Targets</th><th>All 42 motifs</th><th>PPAR-centric (9)</th></tr></thead>
    <tbody>
{rec}    </tbody>
  </table></div>
  <p class="note">With 136 targets, a gain below about 0.04 AUROC cannot be resolved. All 42 motifs sit at chance for PPARD targets (0.506). The PPAR-centric set lifts PPARD and PPARA, within wide intervals.</p>
  <h3>3. Which promoter window carries the signal?</h3>
  <div class="tbl"><table>
    <thead><tr><th>Window around the TSS</th><th>AUROC</th><th>Matched negatives</th></tr></thead>
    <tbody>
{w}    </tbody>
  </table></div>
  <p class="note">The signal sits close to the TSS. A wider window dilutes it. Verified targets are often regulated through PPREs outside the promoter, which a promoter window cannot see.</p>
  <p><strong>Reading:</strong> cutting or weighting motifs does make the hits themselves more reliable. The step that turns hits into a gene score is the weak link. No set met the rule, because no gene-level difference excludes zero.</p>
</section>

'''
p = 'PPRE_motif_prioritisation.html'; s = open(p, encoding='utf8').read()
m = '<section>\n  <h2>The seven ways and what they gave</h2>'
assert m in s; s = s.replace(m, sec + m.replace('The seven ways and what they gave', 'Earlier tests: response and literature genes'), 1)
a = s.index('<section class="verdict">'); b = s.index('</section>', a) + len('</section>')
verdict = '''<section class="verdict">
  <h2>Short answer</h2>
  <ul>
    <li><strong>Yes, cutting or weighting motifs makes the hits more reliable.</strong> Kept hits fall in ChIP peaks clearly more often than removed ones, at equal GC, on datasets not used for the choice. Even sets chosen without outcome data show it.</li>
    <li><strong>That does not turn into a measurable gain for gene-level prediction.</strong> Against 136 verified human PPAR targets, no set differs from all 42 motifs beyond noise. All 42 reach AUROC 0.55, and are at chance for PPARD targets.</li>
    <li><strong>No version met the pre-set rule</strong>, so the score stays on all 42 motifs. The strongest claim is at hit level (site prioritisation). Gene-level scoring is the weak part.</li>
  </ul>
</section>'''
s = s[:a] + verdict + s[b:]
s = s.replace('Seven ways to use fewer or weighted motifs, tested against held-out binding data, hepatocyte response and the published PPRE genes.', 'Seven ways to use fewer or weighted motifs, tested at hit level, against verified PPAR targets, in hepatocytes and on the published PPRE genes.')
s = s.replace('Seven ways to use fewer motifs or weight them, tested against held-out binding data, hepatocyte response and the published PPRE genes.', 'Seven ways to use fewer or weighted motifs, tested at hit level, against verified PPAR targets, in hepatocytes and on the published PPRE genes.')
old = '<li>The score uses all 42 motifs, and all hits are released.</li>'
new = '<li>The score uses all 42 motifs, and all hits are released. Each hit can carry a flag for the reliable subset, since reliable-subset hits are clearly more often bound.</li>'
assert old in s; s = s.replace(old, new, 1)
s = s.replace('(Addenda 5, 10, 11 to 14)', '(Addenda 5, 10, 11 to 15)')
open(p, 'w', encoding='utf8').write(s)
p2 = 'PPRE_v3_to_trainingfree.html'; t = open(p2, encoding='utf8').read()
old2 = 'We tested seven ways to use fewer or weighted motifs. All improve PPARG binding on held-out datasets, none gives a clear gain for gene response, and the cardiomyocyte data cannot separate them from random subsets of the same size. None is adopted.'
new2 = 'We tested seven ways to use fewer or weighted motifs. Kept hits fall in ChIP peaks clearly more often than removed ones, but against 136 verified human PPAR targets no set beats all 42 motifs beyond noise, and the cardiomyocyte data cannot separate them from random subsets of the same size. None is adopted.'
assert old2 in t; t = t.replace(old2, new2, 1); open(p2, 'w', encoding='utf8').write(t)
print("pages updated")
