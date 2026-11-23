"""Locus-number gradient: write every locus of each full orthology pool (results/pool/{set}/) and
subsample lists results/grad/{set}-n{N}r{k}.txt (N in sizes, k = 0..R-1) plus {set}-all.txt.
n50: reuses the 5 disjoint 50-locus subsamples r0..r4 of the replicate experiment (results/loci/{set}-r{k}).
N >= 100: independent draws (overlap allowed) with random.Random(seed + 10000 * N + k).
Usage: select_grad.py <proteomes> <of_out> <busco_out>"""
import sys, os, glob, random, csv, yaml
cfg = yaml.safe_load(open("workflow/config.yaml")); GRAD = yaml.safe_load(open("workflow/config_grad.yaml"))
N, seed, min_occ = cfg["n_loci"], cfg["seed"], cfg["min_occupancy_relaxed"]
prot_dir, of_dir, busco_dir = sys.argv[1:4]
def read_fasta(p):
    d, k = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            k = line[1:].split()[0]; d[k] = []
        elif k: d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}

seqs = {}
for f in glob.glob(f"{prot_dir}/*.fa"):
    seqs.update(read_fasta(f))
taxa = sorted(os.path.basename(f)[:-3] for f in glob.glob(f"{prot_dir}/*.fa"))

pools = {"of_strict": {}, "of_relaxed": {}, "busco": {}}
og_tsv = (glob.glob(f"{of_dir}/**/Orthogroups/Orthogroups.tsv", recursive=True) + glob.glob(f"{of_dir}/**/Orthogroups/MCL_Orthogroups.tsv", recursive=True))[0]
with open(og_tsv) as fh:
    r = csv.reader(fh, delimiter="\t"); header = next(r)
    for row in r:
        og, genes, multi = row[0], {}, False
        for sp, cell in zip(header[1:], row[1:]):
            g = [x.strip() for x in cell.split(",") if x.strip()]
            if len(g) > 1: multi = True
            if len(g) == 1: genes[g[0].split("|")[0]] = g[0]
        if multi: continue
        if len(genes) == len(taxa): pools["of_strict"][og] = genes
        if len(genes) >= min_occ: pools["of_relaxed"][og] = genes
bus = {}
for t in taxa:
    ft = glob.glob(f"{busco_dir}/{t}/run_*/full_table.tsv")[0]
    for line in open(ft):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) > 2 and p[1] == "Complete":
            bus.setdefault(p[0], {})[t] = p[2]
pools["busco"] = {b: g for b, g in bus.items() if len(g) >= min_occ}

os.makedirs("results/grad", exist_ok=True)
with open("results/grad/subsamples.tsv", "w") as out:
    out.write("subsample\tset\tN\treplicate\tpool_size\n")
    for name, pool in pools.items():
        d = f"results/pool/{name}"; os.makedirs(d, exist_ok=True)
        for l, genes in pool.items():
            p = f"{d}/{l}.fa"
            if not os.path.exists(p):
                with open(p, "w") as f:
                    for t, gid in sorted(genes.items()): f.write(f">{t}\n{seqs[gid].replace('*','').replace('U','X')}\n")
        keys = sorted(pool)
        subs = {}
        for k in range(GRAD["n_reps"]):
            for n in GRAD["sizes"]:
                if n == 50:
                    pick = sorted(os.path.basename(x)[:-3] for x in glob.glob(f"results/loci/{name}-r{k}/*.fa"))
                    assert len(pick) == 50 and set(pick) <= set(keys), (name, k)
                else:
                    pick = sorted(random.Random(seed + 10000 * n + k).sample(keys, n))
                subs[(n, k)] = pick
        subs[("all", 0)] = keys
        for (n, k), pick in subs.items():
            sid = f"{name}-all" if n == "all" else f"{name}-n{n}r{k}"
            new = "\n".join(pick) + "\n"
            p = f"results/grad/{sid}.txt"
            if not (os.path.exists(p) and open(p).read() == new): open(p, "w").write(new)
            out.write(f"{sid}\t{name}\t{len(pick) if n == 'all' else n}\t{k}\t{len(pool)}\n")
print(open("results/grad/subsamples.tsv").read().count("\n") - 1, "subsamples")
