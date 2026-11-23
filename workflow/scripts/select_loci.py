"""Step 0: build the three orthology-based locus pools and subsample N loci from each.
of_strict : OrthoFinder orthogroups with exactly one copy in every taxon
of_relaxed: OrthoFinder orthogroups single-copy where present, occupancy >= min_occ
busco     : BUSCO complete single-copy genes, occupancy >= min_occ
Writes results/loci/{set}/{locus}.fa (headers = taxon code) and results/loci/pools.tsv"""
import sys, os, glob, random, csv, yaml
cfg = yaml.safe_load(open("workflow/config.yaml"))
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
    sp_cols = header[1:]
    for row in r:
        og, genes = row[0], {}
        multi = False
        for sp, cell in zip(sp_cols, row[1:]):
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

os.makedirs("results/loci", exist_ok=True)
with open("results/loci/pools.tsv", "w") as out:
    out.write("set\tpool_size\tsampled\tmean_occupancy_sampled\n")
    for name, pool in pools.items():
        rng = random.Random(seed)
        keys = sorted(pool)
        pick = sorted(rng.sample(keys, min(N, len(keys))))
        d = f"results/loci/{name}"; os.makedirs(d, exist_ok=True)
        for k in pick:
            with open(f"{d}/{k}.fa", "w") as f:
                for t, g in sorted(pool[k].items()):
                    f.write(f">{t}\n{seqs[g].replace('*','')}\n")
        occ = sum(len(pool[k]) for k in pick) / max(1, len(pick))
        out.write(f"{name}\t{len(pool)}\t{len(pick)}\t{occ:.1f}\n")
print(open("results/loci/pools.tsv").read())
