#!/usr/bin/env python3
"""Real, column-permuted (5 sets) and reversed copies of the 42 human PPAR/RXR motifs (DESIGN_V4.md Addendum 17).
usage: wsl_make_ctrl42.py <batch_dir> <outdir>"""
import sys, os, re, glob, random
bdir, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
head = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"
mots = []
for f in sorted(glob.glob(os.path.join(bdir, "RXR_PPAR_human42_files_*.meme"))):
    for m in re.split(r"\n(?=MOTIF )", open(f).read())[1:]:
        L = m.rstrip("\n").split("\n"); mots.append((L[0].split()[1], L[1], [l for l in L[2:] if re.match(r"^\s*[0-9.]", l)]))
print(len(mots), "motifs")
def write(name, fn):
    with open(os.path.join(out, name + ".meme"), "w") as fh:
        fh.write(head)
        for mid, hdr, rows in mots:
            suffix, rr = fn(mid, rows); fh.write("\nMOTIF %s%s\n%s\n%s\n" % (mid, suffix, hdr, "\n".join(rr)))
write("real", lambda mid, rows: ("", rows))
for k in range(1, 6):
    rng = random.Random(20261010 * 10 + k)
    def perm(mid, rows, rng=rng, k=k):
        idx = list(range(len(rows)))
        for _ in range(100):
            rng.shuffle(idx)
            if idx != sorted(idx): break
        return "__P%02d" % k, [rows[i] for i in idx]
    write("perm%02d" % k, perm)
write("rev", lambda mid, rows: ("__REV", rows[::-1]))
print("written to", out)
