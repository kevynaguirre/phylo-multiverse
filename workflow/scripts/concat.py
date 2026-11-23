"""Concatenate per-locus alignments (headers = taxon code) into a supermatrix."""
import sys
out, files = sys.argv[1], sys.argv[2:]
def read(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split("|")[0].split()[0]; d[k] = []
        elif k: d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}
alns = [read(f) for f in files]
taxa = sorted(set().union(*[a.keys() for a in alns]))
seqs = {t: [] for t in taxa}
for a in alns:
    L = len(next(iter(a.values())))
    if L == 0: continue
    for t in taxa: seqs[t].append(a.get(t, "-" * L))
with open(out, "w") as f:
    for t in taxa: f.write(f">{t}\n{''.join(seqs[t])}\n")
