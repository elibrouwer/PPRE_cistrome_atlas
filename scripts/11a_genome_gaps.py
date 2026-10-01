"""Chromosome sizes and N-gap intervals of the scanned GRCh38 FASTA (chr1-22, X), used to define where random control sites may fall.

Run in WSL:  python3 11a_genome_gaps.py ~/PPRE_FIMO/Homo_sapiens.GRCh38.dna.primary_assembly.fa out_dir
Writes hg38_chrom_sizes.tsv and hg38_N_gaps.bed (0-based, half-open).
"""
import re, sys, os

fa, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
keep = [str(i) for i in range(1, 23)] + ["X"]
sizes, gaps = {}, []
name, pos, in_n, n_start = None, 0, False, 0


def close_gap(chrom, end):
    global in_n
    if in_n:
        gaps.append((chrom, n_start, end))
        in_n = False


with open(fa, "rb") as fh:
    for line in fh:
        if line.startswith(b">"):
            if name in keep:
                close_gap("chr" + name, pos); sizes["chr" + name] = pos
            name = line[1:].split()[0].decode()
            pos, in_n = 0, False
            continue
        if name not in keep:
            continue
        s = line.strip().upper()
        # N runs inside this line
        for m in re.finditer(rb"N+", s):
            a, b = pos + m.start(), pos + m.end()
            if in_n and m.start() == 0:
                continue                      # gap continues from the previous line
            if in_n and m.start() != 0:
                close_gap("chr" + name, pos)
            if not in_n:
                in_n, n_start = True, a
            if m.end() < len(s):              # gap ends inside this line
                close_gap("chr" + name, b)
        if in_n and not s.endswith(b"N"):
            close_gap("chr" + name, pos)
        pos += len(s)
if name in keep:
    close_gap("chr" + name, pos); sizes["chr" + name] = pos

with open(os.path.join(out, "hg38_chrom_sizes.tsv"), "w") as f:
    for c in ["chr" + k for k in keep]:
        f.write(f"{c}\t{sizes[c]}\n")
with open(os.path.join(out, "hg38_N_gaps.bed"), "w") as f:
    for c, a, b in gaps:
        f.write(f"{c}\t{a}\t{b}\n")
print("chromosomes:", len(sizes), "total bp:", sum(sizes.values()), "| N bp:", sum(b - a for _, a, b in gaps), "| gaps:", len(gaps))
