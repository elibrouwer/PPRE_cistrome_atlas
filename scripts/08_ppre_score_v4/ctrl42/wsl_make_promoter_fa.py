#!/usr/bin/env python3
"""FASTA of the promoter union (cidx, start0, end TSV) from the unmasked GRCh38 primary assembly; names are '<cidx>|<start0>'.
usage: wsl_make_promoter_fa.py <genome.fa> <union.tsv> <out.fa>"""
import sys
ENS = {**{str(i): i - 1 for i in range(1, 23)}, "X": 22, "Y": 23}; INV = {v: k for k, v in ENS.items()}
fa, tsv, out = sys.argv[1:4]
iv = {}
for l in open(tsv):
    c, s, e = l.split()
    if c.isdigit(): iv.setdefault(int(c), []).append((int(s), int(e)))
n = 0
with open(out, "w") as o:
    name = None; buf = []
    def flush():
        global n
        if name in ENS:
            seq = "".join(buf)
            for s, e in sorted(iv.get(ENS[name], [])):
                o.write(">%d|%d\n%s\n" % (ENS[name], s, seq[s:e])); n += 1
    for line in open(fa):
        if line.startswith(">"):
            flush(); name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    flush()
print(n, "windows written")
