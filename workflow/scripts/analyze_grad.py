"""Analyse the locus-number gradient (Snakefile_grad). Outputs in results/analysis_grad/:
  trees.tsv                 one row per tree
  variance_by_N.tsv         per N: R2 orthology, subsample-within-orthology, inference (partial), subsample x inference,
                            residual; nested permutation test of orthology against subsamples
  rf_by_N.tsv               per N and method: mean RF between replicate subsamples of the same set (sampling noise),
                            mean RF to the full-pool tree of the same set and method, mean RF between sets
  full_pool_rf.tsv          RF among the 12 full-pool trees (3 orthology sets x 4 methods)
  clades_by_N.tsv           frequency of each variable clade per N x orthology set (+ presence in full-pool trees)
  fig_grad_variance.png, fig_grad_rf.png
"""
import glob, os, itertools, yaml
import numpy as np, pandas as pd, dendropy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.ticker

cfg = yaml.safe_load(open("workflow/config.yaml")); OUT = cfg["outgroup"]
OD = os.environ.get("GRAD_OUT", "results/analysis_grad"); os.makedirs(OD, exist_ok=True)
meta = pd.read_csv("results/grad/subsamples.tsv", sep="\t").set_index("subsample")
rows = []
for f in sorted(glob.glob("results/gctree/*/*.treefile")):
    s, m = f.split("/")[2], os.path.basename(f)[:-9]; rows.append(dict(sub=s, inference=f"concat_{m}", file=f))
for f in sorted(glob.glob("results/gastral/*/*.tre")):
    s, m = f.split("/")[2], os.path.basename(f)[:-4]; rows.append(dict(sub=s, inference=f"astral_{m}", file=f))
df = pd.DataFrame(rows)
df["orthology"] = df["sub"].map(meta.set); df["full"] = df["sub"].str.endswith("-all")
df["N"] = np.where(df.full, "all", df["sub"].map(meta.N).astype(str)); df["replicate"] = df["sub"].map(meta.replicate)
tns = dendropy.TaxonNamespace(); trees, keep = [], []
for i, r in df.iterrows():
    s = open(r.file).read().strip()
    if not s: continue
    trees.append(dendropy.Tree.get(data=s, schema="newick", taxon_namespace=tns, preserve_underscores=True)); keep.append(i)
df = df.loc[keep].reset_index(drop=True)
taxa = sorted(t.label for t in tns); ntax = len(taxa); maxrf = 2 * (ntax - 3)
def clades(tr):
    tr.encode_bipartitions(); out = set()
    for b in tr.bipartition_encoding:
        side = frozenset(t.label for t in b.leafset_taxa(tns))
        if OUT in side: side = frozenset(taxa) - side
        if 1 < len(side) < ntax - 1: out.add(side)
    return out
C = [clades(t) for t in trees]; n = len(C)
D = np.zeros((n, n))
for i, j in itertools.combinations(range(n), 2): D[i, j] = D[j, i] = len(C[i] ^ C[j]) / maxrf
df.drop(columns="file").to_csv(f"{OD}/trees.tsv", sep="\t", index=False)

def gower(Dm):
    m = len(Dm); J = np.eye(m) - 1 / m; return J @ (-0.5 * Dm ** 2) @ J
def ss(G, d, fs):
    X = [np.ones(len(d))]
    for f in fs:
        for l in sorted(d[f].unique())[1:]: X.append((d[f] == l).astype(float).values)
    X = np.column_stack(X); H = X @ np.linalg.pinv(X); return np.trace(H @ G @ H)
rng = np.random.default_rng(1)
sizes = sorted(int(x) for x in df.N.unique() if x != "all")
vrows, rfrows = [], []
full_tree = {(r.orthology, r.inference): i for i, r in df[df.full].iterrows()}
for N in sizes:
    ix = np.where(df.N == str(N))[0]; d = df.loc[ix].reset_index(drop=True); G = gower(D[np.ix_(ix, ix)]); SST = np.trace(G)
    d["sxi"] = d["sub"] + ":" + d.inference
    so, ssub = ss(G, d, ["orthology"]), ss(G, d, ["sub"]); sfull = ss(G, d, ["sub", "inference"])
    a, R = d.orthology.nunique(), d.replicate.nunique()
    uo = d.groupby("sub").orthology.first()
    def F(lbl):
        s_ = ss(G, d.assign(orthology=d["sub"].map(lbl)), ["orthology"]); return (s_ / (a - 1)) / ((ssub - s_) / (a * (R - 1)))
    obs = F(uo); P = 4999
    p = (1 + sum(F(pd.Series(rng.permutation(uo.values), index=uo.index)) >= obs for _ in range(P))) / (P + 1)
    vrows.append(dict(N=N, n_trees=len(d), R2_orthology=so / SST, p_orthology_nested=p, R2_subsample_within_orth=(ssub - so) / SST,
                      R2_inference_partial=(sfull - ssub) / SST, R2_inference_marginal=ss(G, d, ["inference"]) / SST,
                      R2_subsample_x_inference=(ss(G, d, ["sxi"]) - sfull) / SST, mean_RF_all_pairs=D[np.ix_(ix, ix)].sum() / (len(ix) * (len(ix) - 1))))
    for inf, g in df.loc[ix].groupby("inference"):
        same = [D[i, j] for i, j in itertools.combinations(g.index, 2) if df.orthology[i] == df.orthology[j]]
        diff = [D[i, j] for i, j in itertools.combinations(g.index, 2) if df.orthology[i] != df.orthology[j]]
        tofull = [D[i, full_tree[(df.orthology[i], inf)]] for i in g.index if (df.orthology[i], inf) in full_tree]
        rfrows.append(dict(N=N, inference=inf, RF_between_replicates_same_set=np.mean(same), RF_between_sets=np.mean(diff),
                           RF_to_full_pool_tree=np.mean(tofull) if tofull else np.nan))
vp = pd.DataFrame(vrows); vp.to_csv(f"{OD}/variance_by_N.tsv", sep="\t", index=False)
rf = pd.DataFrame(rfrows); rf.to_csv(f"{OD}/rf_by_N.tsv", sep="\t", index=False)
fx = df[df.full].index
fp = pd.DataFrame(D[np.ix_(fx, fx)], index=df.orthology[fx] + "|" + df.inference[fx], columns=df.orthology[fx] + "|" + df.inference[fx])
fp.to_csv(f"{OD}/full_pool_rf.tsv", sep="\t")

allc = {}
for i, cs in enumerate(C):
    for c in cs: allc.setdefault(c, set()).add(i)
rec = []
for c, ix in allc.items():
    pres = np.array([1.0 if i in ix else 0.0 for i in range(n)])
    if not (0.03 < pres.mean() < 0.97): continue
    d = dict(clade=",".join(sorted(c)), freq_all_trees=pres.mean())
    for (N, o), g in df.assign(p=pres).groupby(["N", "orthology"]): d[f"N{N}|{o}"] = g.p.mean()
    rec.append(d)
cl = pd.DataFrame(rec).sort_values("freq_all_trees", ascending=False)
cols = ["clade", "freq_all_trees"] + sorted([c for c in cl.columns if c.startswith("N") and "|" in c and c != "clade"],
        key=lambda c: (c.split("|")[1], 10**9 if c.split("|")[0][1:] == "all" else int(c.split("|")[0][1:])))
cl = cl[cols]; cl.to_csv(f"{OD}/clades_by_N.tsv", sep="\t", index=False)

fig, ax = plt.subplots(figsize=(6, 3.6))
for col, lab, c in [("R2_subsample_within_orth", "subsample (within orthology)", "#72B7B2"), ("R2_subsample_x_inference", "subsample × method", "#E45756"),
                    ("R2_inference_partial", "method (main effect)", "#F58518"), ("R2_orthology", "orthology", "#4C78A8")]:
    ax.plot(vp.N, vp[col].clip(lower=0), "o-", color=c, label=lab)
ax.set_xscale("log"); ax.set_xticks(sizes, [str(s) for s in sizes]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlabel("loci per subsample"); ax.set_ylabel("share of topological variance")
ax.set_ylim(0, 1); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(f"{OD}/fig_grad_variance.png", dpi=200)
fig, ax = plt.subplots(figsize=(6, 3.6)); r_ = rf.groupby("N")[["RF_between_replicates_same_set", "RF_between_sets", "RF_to_full_pool_tree"]].mean()
for col, lab, c in [("RF_between_replicates_same_set", "between subsamples, same set", "#72B7B2"), ("RF_between_sets", "between orthology sets", "#4C78A8"),
                    ("RF_to_full_pool_tree", "to full-pool tree (same set, method)", "#333333")]:
    ax.plot(r_.index, r_[col], "o-", color=c, label=lab)
ax.set_xscale("log"); ax.set_xticks(sizes, [str(s) for s in sizes]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlabel("loci per subsample"); ax.set_ylabel("mean normalised RF")
ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(f"{OD}/fig_grad_rf.png", dpi=200)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
print(vp.to_string(index=False)); print(rf.to_string(index=False)); print(fp.round(3).to_string())
print(cl.to_string(index=False))
