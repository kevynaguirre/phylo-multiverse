"""Normalise trimmer output: keep original taxon names/order, drop trailing annotations."""
import sys
orig, tmp, out = sys.argv[1:4]
def read(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]; d[k] = []
        elif k: d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}
o, t = read(orig), read(tmp)
L = len(next(iter(t.values()))) if t else 0
with open(out, "w") as f:
    for k in o:
        f.write(f">{k}\n{t.get(k, '-' * L)}\n")
