"""Step 0b: replicate locus subsamples per orthology set.
Builds the same three pools as select_loci.py, then draws R+1 DISJOINT subsamples of N loci per set:
  r0 = exactly the subsample of the original pilot (random.Random(seed).sample, seed from config.yaml)
  r1..rR = successive non-overlapping draws from the remaining loci (random.Random(seed + 1000 * k))
Writes results/loci/{set}-r{k}/{locus}.fa and results/loci/pools_reps.tsv
Usage: select_loci_reps.py <proteomes> <of_out> <busco_out> <R>"""
import sys, os, glob, random, csv, yaml
cfg = yaml.safe_load(open("workflow/config.yaml"))
N, seed, min_occ = cfg["n_loci"], cfg["seed"], cfg["min_occupancy_relaxed"]
prot_dir, of_dir, busco_dir, R = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

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

os.makedirs("results/loci", exist_ok=True)
with open("results/loci/pools_reps.tsv", "w") as out:
    out.write("set\treplicate\tpool_size\tsampled\tmean_occupancy_sampled\n")
    for name, pool in pools.items():
        keys = sorted(pool)
        draws = [sorted(random.Random(seed).sample(keys, min(N, len(keys))))]
        used = set(draws[0])
        for k in range(1, R + 1):
            rest = sorted(set(keys) - used)
            pick = sorted(random.Random(seed + 1000 * k).sample(rest, min(N, len(rest))))
            draws.append(pick); used |= set(pick)
        for k, pick in enumerate(draws):
            d = f"results/loci/{name}-r{k}"; os.makedirs(d, exist_ok=True)
            for l in pick:
                with open(f"{d}/{l}.fa", "w") as f:
                    for t, g in sorted(pool[l].items()):
                        f.write(f">{t}\n{seqs[g].replace('*','')}\n")
            occ = sum(len(pool[l]) for l in pick) / max(1, len(pick))
            out.write(f"{name}\tr{k}\t{len(pool)}\t{len(pick)}\t{occ:.1f}\n")
print(open("results/loci/pools_reps.tsv").read())
