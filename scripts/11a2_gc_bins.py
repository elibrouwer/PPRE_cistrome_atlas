"""GC and N counts per 1 kb window of the scanned GRCh38 FASTA (chr1-22, X), used for the GC-matched control in 11b.

Run in WSL:  python3 11a2_gc_bins.py ~/PPRE_FIMO/Homo_sapiens.GRCh38.dna.primary_assembly.fa out_dir
Writes gc_bins_1kb.tsv with columns chrom, bin (0-based, window = bin*1000 .. bin*1000+999), gc, n, length.
"""
import sys, os

fa, out = sys.argv[1], sys.argv[2]
keep = [str(i) for i in range(1, 23)] + ["X"]
W = 1000
rows = 0
with open(fa, "rb") as fh, open(os.path.join(out, "gc_bins_1kb.tsv"), "w") as fo:
    fo.write("chrom\tbin\tgc\tn\tlength\n")
    name, buf = None, []

    def flush():
        global rows
        if name in keep:
            seq = b"".join(buf).upper()
            for b in range(0, len(seq), W):
                s = seq[b:b + W]
                fo.write(f"chr{name}\t{b // W}\t{s.count(b'G') + s.count(b'C')}\t{s.count(b'N')}\t{len(s)}\n")
                rows += 1

    for line in fh:
        if line.startswith(b">"):
            flush()
            name, buf = line[1:].split()[0].decode(), []
        elif name in keep:
            buf.append(line.strip())
    flush()
print("windows written:", rows)
