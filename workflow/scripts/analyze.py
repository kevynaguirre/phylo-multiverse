"""Analyse the multiverse of species trees.
Outputs (results/analysis/):
  trees.tsv            one row per pipeline (factors + RF to majority consensus)
  rf_matrix.tsv        pairwise normalised RF distances
  variance_partition.tsv  marginal and partial R2 of each factor (distance-based, Gower/dbRDA)
  clades.tsv           every clade seen in >=5% of trees: frequency (=1-DSI), per-factor-level frequency,
                       per-factor effect (max-min freq across levels), mean gCF/sCF when present
  figures: fig_variance.png, fig_clade_heatmap.png, fig_support_vs_stability.png
"""
import glob, os, sys, itertools, re
import numpy as np, pandas as pd, dendropy, yaml
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

cfg = yaml.safe_load(open("workflow/config.yaml"))
OUT = cfg["outgroup"]; os.makedirs("results/analysis", exist_ok=True)
tns = dendropy.TaxonNamespace()
rows, trees = [], []
for f in sorted(glob.glob("results/ctree/*/*/*/*.treefile")):
    o, a, t, m = f.split("/")[2:6]; m = m.replace(".treefile", "")
    rows.append(dict(orthology=o, aligner=a, trimming=t, method="concat", model=m, file=f))
for f in sorted(glob.glob("results/astral/*/*/*/*.tre")):
    o, a, t, m = f.split("/")[2:6]; m = m.replace(".tre", "")
    rows.append(dict(orthology=o, aligner=a, trimming=t, method="astral", model=m, file=f))
df = pd.DataFrame(rows)
df["inference"] = df.method + "_" + df.model
keep = []
for i, r in df.iterrows():
    s = open(r.file).read().strip()
    if not s: continue
    tr = dendropy.Tree.get(data=s, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    trees.append(tr); keep.append(i)
df = df.loc[keep].reset_index(drop=True)
taxa = sorted(t.label for t in tns)
ntax = len(taxa)

def clades(tr):
    """Non-trivial clades, expressed on the side NOT containing the outgroup (rooting-independent)."""
    tr.encode_bipartitions()
    out = set()
    for b in tr.bipartition_encoding:
        side = frozenset(t.label for t in b.leafset_taxa(tns))
        if OUT in side: side = frozenset(taxa) - side
        if 1 < len(side) < ntax - 1: out.add(side)
    return out
C = [clades(t) for t in trees]
n = len(C)
# pairwise normalised RF (unrooted)
D = np.zeros((n, n))
maxrf = 2 * (ntax - 3)
for i, j in itertools.combinations(range(n), 2):
    D[i, j] = D[j, i] = len(C[i] ^ C[j]) / maxrf
pd.DataFrame(D).to_csv("results/analysis/rf_matrix.tsv", sep="\t")

# --- distance-based variance partitioning (Gower-centred, dbRDA; McArdle & Anderson 2001)
A = -0.5 * D ** 2
J = np.eye(n) - 1 / n
G = J @ A @ J
SST = np.trace(G)
FACT = ["orthology", "aligner", "trimming", "inference"]
def design(fs):
    X = [np.ones(n)]
    for f in fs:
        lv = sorted(df[f].unique())
        for l in lv[1:]: X.append((df[f] == l).astype(float).values)
    return np.column_stack(X)
def ss(fs):
    X = design(fs); H = X @ np.linalg.pinv(X)
    return np.trace(H @ G @ H)
rng = np.random.default_rng(1)
vp = []
full = ss(FACT)
for f in FACT:
    marg = ss([f]) / SST
    others = [g for g in FACT if g != f]
    part = (full - ss(others)) / SST
    # permutation test of the partial effect: permute factor labels
    obs = part; cnt = 0; P = 999
    for _ in range(P):
        perm = df[f].sample(frac=1, random_state=int(rng.integers(1e9))).values
        old = df[f].copy(); df[f] = perm
        cnt += ((ss(FACT) - ss(others)) / SST) >= obs
        df[f] = old
    vp.append(dict(factor=f, marginal_R2=marg, partial_R2=part, p_perm=(cnt + 1) / (P + 1)))
vp.append(dict(factor="all_factors", marginal_R2=full / SST, partial_R2=np.nan, p_perm=np.nan))
vp = pd.DataFrame(vp); vp.to_csv("results/analysis/variance_partition.tsv", sep="\t", index=False)

# --- clade-level stability
allc = {}
for i, cs in enumerate(C):
    for c in cs: allc.setdefault(c, []).append(i)
majority = {c for c, ix in allc.items() if len(ix) / n > 0.5}
df["RF_to_majority"] = [len(cs ^ majority) / maxrf for cs in C]
df["n_trees"] = n
df.drop(columns="file").to_csv("results/analysis/trees.tsv", sep="\t", index=False)

# concordance factors from IQ-TREE (concat trees): map clade -> gCF/sCF lists
cf = {}
for i, r in df[df.method == "concat"].iterrows():
    stat = r.file.replace(".treefile", ".cf.cf.stat"); tf = r.file.replace(".treefile", ".cf.cf.tree")
    if not (os.path.exists(stat) and os.path.exists(tf)): continue
    tr = dendropy.Tree.get(path=tf, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    # IQ-TREE writes concordance factors as internal-node labels "gCF/sCF"
    for nd in tr.postorder_internal_node_iter():
        if nd.label is None or nd is tr.seed_node: continue
        parts = str(nd.label).split("/")
        if len(parts) < 2: continue
        try: g, s_ = float(parts[0]), float(parts[1])
        except ValueError: continue
        side = frozenset(l.taxon.label for l in nd.leaf_iter())
        if OUT in side: side = frozenset(taxa) - side
        if not (1 < len(side) < ntax - 1): continue
        cf.setdefault(side, []).append((g, s_))

# ASTRAL local posterior probabilities (node labels of ASTRAL trees)
lpp = {}
for i, r in df[df.method == "astral"].iterrows():
    tr = dendropy.Tree.get(path=r.file, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    for nd in tr.postorder_internal_node_iter():
        if nd.label is None or nd is tr.seed_node: continue
        try: v = float(nd.label)
        except ValueError: continue
        side = frozenset(l.taxon.label for l in nd.leaf_iter())
        if OUT in side: side = frozenset(taxa) - side
        if 1 < len(side) < ntax - 1: lpp.setdefault(side, []).append(v)

rec = []
for c, ix in allc.items():
    freq = len(ix) / n
    if freq < 0.05: continue
    d = dict(clade=",".join(sorted(c)), size=len(c), freq=freq, DSI=1 - freq, in_majority=c in majority)
    pres = np.zeros(n); pres[ix] = 1
    for f in FACT:
        lv = df.groupby(f).apply(lambda g: pres[g.index].mean())
        for l, v in lv.items(): d[f"{f}={l}"] = v
        d[f"effect_{f}"] = lv.max() - lv.min()
    if c in cf:
        d["mean_gCF"] = np.mean([x[0] for x in cf[c]]); d["mean_sCF"] = np.mean([x[1] for x in cf[c]])
    if c in lpp: d["mean_localPP"] = np.mean(lpp[c])
    rec.append(d)
cl = pd.DataFrame(rec).sort_values(["in_majority", "freq"], ascending=[False, True])
cl.to_csv("results/analysis/clades.tsv", sep="\t", index=False)

# --- figures
fig, ax = plt.subplots(figsize=(6, 3.5))
v = vp[vp.factor != "all_factors"]
x = np.arange(len(v)); ax.bar(x - 0.2, v.marginal_R2, 0.4, label="marginal R²", color="#4C78A8")
ax.bar(x + 0.2, v.partial_R2, 0.4, label="partial R²", color="#F58518")
ax.set_xticks(x, v.factor); ax.set_ylabel("Share of topological variance (RF)")
ax.legend(frameon=False); ax.set_title(f"{n} species trees, {ntax} taxa"); fig.tight_layout()
fig.savefig("results/analysis/fig_variance.png", dpi=200)

var = cl[(cl.freq < 0.98)].copy()
if len(var):
    order = df.sort_values(FACT).index
    M = np.array([[1 if (frozenset(r.split(",")) in C[i]) else 0 for i in order] for r in var.clade])
    fig, ax = plt.subplots(figsize=(12, max(2.5, 0.3 * len(var) + 1.5)))
    ax.imshow(M, aspect="auto", cmap="Greys", interpolation="nearest")
    lab = [f"{r.clade[:60]} ({r.freq:.2f})" for r in var.itertuples()]
    ax.set_yticks(range(len(var)), lab, fontsize=6)
    ax.set_xlabel("pipelines (sorted by orthology > aligner > trimming > inference)")
    fig.tight_layout(); fig.savefig("results/analysis/fig_clade_heatmap.png", dpi=200)

if "mean_gCF" in cl:
    s = cl.dropna(subset=["mean_gCF"])
    cols = [c_ for c_ in ["mean_gCF", "mean_sCF", "mean_localPP"] if c_ in s]
    fig, ax = plt.subplots(1, len(cols), figsize=(4 * len(cols), 3.5))
    for a_, col in zip(ax, cols):
        a_.scatter(s[col], s.freq, s=14, color="#4C78A8"); a_.set_xlabel(col); a_.set_ylabel("frequency across multiverse (1-DSI)")
    fig.tight_layout(); fig.savefig("results/analysis/fig_support_vs_stability.png", dpi=200)
print(vp.to_string(index=False))
print(cl[["clade", "freq", "in_majority"] + [c for c in cl.columns if c.startswith("effect_")]].head(30).to_string(index=False))
