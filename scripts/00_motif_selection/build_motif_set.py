"""Build the complete PPAR / RXR / PPAR::RXR motif set (JASPAR 2026, HOCOMOCO CORE, CIS-BP 3.10 direct, human)."""
import csv, io, json, os, re, shutil, tarfile, time, urllib.request, urllib.parse

SP = os.path.dirname(os.path.abspath(__file__))
ROOT = r"C:\Users\brouw\OneDrive - Universiteit Utrecht\Master Bioinformatics\Minor internship\PPAR_RXR_motifs_Hocomoco_Jaspar"
OUT = os.path.join(ROOT, "Motifs_update_2026-09-30")
for d in ("JASPAR", "HOCOMOCO", "CIS-BP"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)

GENES = ["PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG"]


def fetch(url, tries=5, binary=False):
    err = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (academic motif retrieval)"})
            with urllib.request.urlopen(req, timeout=180) as r:
                b = r.read()
                return b if binary else b.decode("utf-8")
        except Exception as e:  # noqa
            err = e
            time.sleep(2 + 2 * i)
    raise RuntimeError(f"{url}: {err}")


def iupac(rows):
    code = {"A": "A", "C": "C", "G": "G", "T": "T", "AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K",
            "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B"}
    s = ""
    for p in rows:
        m = max(p)
        if m >= 0.5:
            s += "ACGT"[p.index(m)]
            continue
        k = "".join(b for b, v in zip("ACGT", p) if v >= 0.25)
        s += code.get(k, "N") if k else "N"
    return s


def meme_block(name, rows, nsites, url=None):
    t = f"MOTIF {name}\nletter-probability matrix: alength= 4 w= {len(rows)} nsites= {int(round(nsites))}\n"
    t += "\n".join(" ".join(f"{v:.6f}" for v in r) for r in rows) + "\n"
    if url:
        t += f"URL {url}\n"
    return t + "\n"


MEME_HEAD = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"

current22 = set("""MA1148.1 MA1148.2 MA0065.1 MA1550.1 MA1550.2 MA0066.1 MA0066.2 MA0855.1 MA1555.1 MA0856.1 MA1556.1
PPARA.H13CORE.0.P.B PPARA.H13CORE.1.P.B PPARD.H13CORE.0.PSM.A PPARG.H13CORE.0.P.B PPARG.H13CORE.1.P.B
RXRA.H13CORE.0.PS.A RXRA.H13CORE.1.S.C RXRA.H13CORE.2.SM.B RXRA.H13CORE.3.P.B RXRB.H13CORE.0.P.C RXRB.H13CORE.1.SM.B""".split())
cisbp_curated8 = set("M12541 M05794 M05796 M12496 M12497 M12498 M03621 M03623".split())

inventory, all_meme = [], {}

# ---------------------------------------------------------------- JASPAR
det = json.load(open(os.path.join(SP, "jaspar_details_raw.json")))
keep = {}
for mid, x in det.items():
    nm = x["name"].upper()
    parts = nm.split("::")
    has_ppar = any(p.startswith("PPAR") for p in parts)
    has_rxr = any(p.startswith("RXR") for p in parts)
    if (has_ppar or has_rxr) and all(p.startswith(("PPAR", "RXR")) for p in parts):
        keep[mid] = x
print("JASPAR kept:", len(keep))
for mid in sorted(keep):
    x = keep[mid]
    txt = fetch(f"https://jaspar.elixir.no/api/v1/matrix/{mid}/?format=jaspar")
    open(os.path.join(OUT, "JASPAR", f"{mid}.jaspar"), "w", newline="\n").write(txt)
    pfm = x["pfm"]
    cols = list(zip(*[pfm[b] for b in "ACGT"]))
    rows = [[v / sum(c) for v in c] for c in cols]
    nsites = sum(cols[0])
    all_meme[mid] = meme_block(mid, rows, nsites, f"https://jaspar.elixir.no/matrix/{mid}/")
    open(os.path.join(OUT, "JASPAR", f"{mid}_meme_format.meme"), "w", newline="\n").write(
        MEME_HEAD + meme_block(mid, rows, nsites, f"https://jaspar.elixir.no/matrix/{mid}/"))
    species = "; ".join(s["name"] for s in x["species"])
    parts = x["name"].upper().split("::")
    cat = "PPAR:RXR heterodimer" if len(parts) == 2 else ("PPAR" if parts[0].startswith("PPAR") else "RXR")
    inventory.append(dict(
        database="JASPAR 2026", name=x["name"], motif_id=mid, native_id=mid, category=cat,
        data_type=(x.get("type") or "").strip() or "ChIP-seq (ReMap/GTRD-derived)" if not x.get("type") else x["type"],
        source_species=species, width=len(rows), consensus=iupac(rows),
        status="in current 22" if mid in current22 else "NEW",
        source=f"PMID {'; '.join(x.get('pubmed_ids') or [])}", notes=(x.get("comment") or "").strip()[:200]))

# ---------------------------------------------------------------- HOCOMOCO (CORE, human bundle; v14 matrices == v13)
hb = "https://hocomoco14.autosome.org/final_bundle/hocomoco13/H13CORE"
ann = {}
for l in open(os.path.join(SP, "hocomoco", "H13CORE_annotation.jsonl"), encoding="utf-8"):
    r = json.loads(l)
    if r["tf"] in GENES:
        ann[r["name"]] = r
print("HOCOMOCO H13CORE PPAR/RXR:", len(ann))
tars = {}
for kind in ("pcm", "pwm"):
    b = fetch(f"{hb}/H13CORE_{kind}.tar.gz", binary=True)
    tf = tarfile.open(fileobj=io.BytesIO(b))
    tars[kind] = {os.path.basename(m.name): tf.extractfile(m).read().decode() for m in tf.getmembers() if m.isfile()}
print("tar members example:", list(tars["pcm"])[:3])
memebundle = fetch(f"{hb}/formatted_motifs/H13CORE_meme_format.meme")
blocks = re.split(r"(?=^MOTIF )", memebundle, flags=re.M)
meme_by = {b.split()[1]: b for b in blocks if b.startswith("MOTIF ")}
LET = {"P": "ChIP-Seq", "S": "HT-SELEX", "M": "Methyl-HT-SELEX", "G": "GHT-SELEX", "B": "PBM", "L": "SMiLE-seq"}
for nm in sorted(ann):
    r = ann[nm]
    for kind in ("pcm", "pwm"):
        key = next(k for k in tars[kind] if k.startswith(nm + "."))
        open(os.path.join(OUT, "HOCOMOCO", f"{nm}.{kind}"), "w", newline="\n").write(tars[kind][key])
    mb = meme_by[nm]
    open(os.path.join(OUT, "HOCOMOCO", f"{nm}_meme_format.meme"), "w", newline="\n").write(MEME_HEAD + mb.strip() + "\n")
    rows = r["pfm"]
    all_meme[nm] = mb.strip() + "\n\n"
    dt = " + ".join(LET.get(c, c) for c in r["datatype"])
    inventory.append(dict(
        database="HOCOMOCO v13/v14 CORE", name=r["tf"], motif_id=nm, native_id=nm, category="PPAR" if r["tf"].startswith("PPAR") else "RXR",
        data_type=dt, source_species=f"{r['original_motif']['species'].title()} (TF: human; v14 matrix identical)",
        width=r["length"], consensus=iupac(rows), status="in current 22" if nm in current22 else "NEW",
        source=r["original_motif"]["origin"], notes=f"quality {r['quality']}; {r['num_words']} sites"))

# ---------------------------------------------------------------- CIS-BP (direct, human, matrices published)
cdir = os.path.join(SP, "cisbp", "human")
rows_ = list(csv.reader(open(os.path.join(cdir, "TF_Information_all_motifs_plus.txt"), encoding="utf-8", errors="replace"), delimiter="\t"))
h = rows_[0]
ix = lambda c: h.index(c)
seen = set()
cis = []
for r in rows_[1:]:
    if r[ix("TF_Name")] in GENES and r[ix("TF_Status")] == "D" and r[ix("Motif_ID")] != ".":
        if r[ix("Motif_ID")] in seen:
            continue
        seen.add(r[ix("Motif_ID")])
        cis.append(r)
skipped = []
for r in sorted(cis, key=lambda r: (r[ix("TF_Name")], r[ix("Motif_ID")])):
    mid, gene = r[ix("Motif_ID")], r[ix("TF_Name")]
    src = open(os.path.join(cdir, "pwms_all_motifs", f"{mid}.txt")).read()
    lines = [l.split() for l in src.strip().splitlines()[1:]]
    if not lines:
        skipped.append((gene, mid, r[ix("MSource_Identifier")]))
        continue
    rowsp = [[float(v) for v in l[1:5]] for l in lines]
    name = f"{gene}_{mid}"
    open(os.path.join(OUT, "CIS-BP", f"{mid}.txt"), "w", newline="\n").write(src)
    all_meme[name] = meme_block(name, rowsp, 20, f"https://cisbp.ccbr.utoronto.ca/TFreport.php?searchTF={mid}")
    open(os.path.join(OUT, "CIS-BP", f"{name}_meme_format.meme"), "w", newline="\n").write(MEME_HEAD + all_meme[name])
    typ = r[ix("MSource_Type")]
    typ = {"Misc": "HOCOMOCO v11 / GHT-SELEX / HT-SELEX (Misc)"}.get(typ, typ)
    inventory.append(dict(
        database="CIS-BP 3.10", name=gene, motif_id=name, native_id=f"{mid} ({r[ix('MSource_Identifier')]})", category="PPAR" if gene.startswith("PPAR") else "RXR",
        data_type=f"{r[ix('MSource_Type')]}", source_species="Human TF (direct evidence); species of source data not recorded by CIS-BP",
        width=len(rowsp), consensus=iupac(rowsp), status="in current 30 (curated CIS-BP)" if mid.split("_")[0] in cisbp_curated8 else "NEW",
        source=f"{r[ix('MSource_Author')]} {r[ix('MSource_Year')]}; PMID {r[ix('PMID')]}",
        notes="nsites placeholder 20 in MEME" + ("; labelled 'dimer' by source" if "dimer" in r[ix("MSource_Identifier")].lower() else "")))
print("CIS-BP written:", len([i for i in inventory if i["database"].startswith("CIS")]), "skipped (no published matrix):", skipped)

# ---------------------------------------------------------------- combined outputs
open(os.path.join(OUT, "PPAR_RXR_all_motifs_combined.meme"), "w", newline="\n").write(MEME_HEAD + "".join(all_meme[k] for k in all_meme))
order = {"JASPAR 2026": 0, "HOCOMOCO v13/v14 CORE": 1, "CIS-BP 3.10": 2}
inventory.sort(key=lambda d: (d["name"].split("::")[0][:6] if False else d["name"].upper(), order[d["database"]], d["motif_id"]))
with open(os.path.join(OUT, "motif_inventory.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(inventory[0].keys()))
    w.writeheader()
    w.writerows(inventory)
print("total motifs:", len(inventory), "| in combined meme:", len(all_meme))
