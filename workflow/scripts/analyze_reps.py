"""Separate the orthology effect from the locus-sampling effect using replicate subsamples.
Input: species trees of the replicate design (config_reps.yaml): {set}-r{k} x aligners x trimmers x 5 inference methods.
Outputs (results/analysis_reps/):
  trees.tsv, rf_matrix.tsv
  variance_partition.tsv   nested distance-based partition (Gower/dbRDA): orthology, subsample within
                           orthology, trimming, inference, residual; permutation tests with restricted permutations
  per_replicate_R2.tsv     R2 of orthology and inference computed inside each replicate (as in the single-subsample pilot)
  null_orthology_R2.tsv    R2 of a fake 'orthology' factor made of 3 subsamples of the SAME orthology set
  rf_within_between.tsv    mean RF between trees built with the same method on different subsamples,
                           same orthology set vs different orthology sets
  clades_by_subsample.tsv  frequency of each variable clade in every set x replicate
  locus_overlap.tsv        overlap of sampled loci (via Scer/any shared protein IDs) between subsamples
  figures fig_reps_*.png
"""
import glob, os, itertools, re, yaml
import numpy as np, pandas as pd, dendropy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

cfg = yaml.safe_load(open("workflow/config.yaml")); cfg.update(yaml.safe_load(open("workflow/config_reps.yaml")))
OUT = cfg["outgroup"]; OD = os.environ.get("REPS_OUT", "results/analysis_reps"); os.makedirs(OD, exist_ok=True)
REPS, ALN, TRIM = cfg["replicates"], cfg["aligners"], cfg["trimmers"]
rows = []
for f in sorted(glob.glob("results/ctree/*-r*/*/*/*.treefile")):
    o, a, t, m = f.split("/")[2:6]; rows.append(dict(sub=o, aligner=a, trimming=t, method="concat", model=m[:-9], file=f))
for f in sorted(glob.glob("results/astral/*-r*/*/*/*.tre")):
    o, a, t, m = f.split("/")[2:6]; rows.append(dict(sub=o, aligner=a, trimming=t, method="astral", model=m[:-4], file=f))
df = pd.DataFrame(rows)
df["orthology"] = df["sub"].str.replace(r"-r\d+$", "", regex=True)
df["replicate"] = df["sub"].str.extract(r"-r(\d+)$")[0].astype(int)
df["inference"] = df.method + "_" + df.model
EXCL = [x for x in os.environ.get("REPS_EXCLUDE", "").split(",") if x]
df = df[~(df.method + "_" + df.model).isin(EXCL)]
df = df[df.replicate.isin(REPS) & df.aligner.isin(ALN) & df.trimming.isin(TRIM)].reset_index(drop=True)
tns = dendropy.TaxonNamespace(); trees, keep = [], []
for i, r in df.iterrows():
    s = open(r.file).read().strip()
    if not s: continue
    trees.append(dendropy.Tree.get(data=s, schema="newick", taxon_namespace=tns, preserve_underscores=True)); keep.append(i)
df = df.loc[keep].reset_index(drop=True)
taxa = sorted(t.label for t in tns); ntax = len(taxa)
cells = df.groupby(["orthology", "replicate", "trimming", "inference"]).size()
print(f"{len(df)} trees; design cells with a tree: {len(cells)}; expected "
      f"{df.orthology.nunique() * len(REPS) * len(TRIM) * df.inference.nunique() * len(ALN)}")

def clades(tr):
    tr.encode_bipartitions(); out = set()
    for b in tr.bipartition_encoding:
        side = frozenset(t.label for t in b.leafset_taxa(tns))
        if OUT in side: side = frozenset(taxa) - side
        if 1 < len(side) < ntax - 1: out.add(side)
    return out
C = [clades(t) for t in trees]; n = len(C); maxrf = 2 * (ntax - 3)
D = np.zeros((n, n))
for i, j in itertools.combinations(range(n), 2):
    D[i, j] = D[j, i] = len(C[i] ^ C[j]) / maxrf
pd.DataFrame(D).to_csv(f"{OD}/rf_matrix.tsv", sep="\t")

def gower(Dm):
    m = len(Dm); J = np.eye(m) - 1 / m
    return J @ (-0.5 * Dm ** 2) @ J
def ss(G, d, fs):
    X = [np.ones(len(d))]
    for f in fs:
        for l in sorted(d[f].unique())[1:]: X.append((d[f] == l).astype(float).values)
    X = np.column_stack(X); H = X @ np.linalg.pinv(X)
    return np.trace(H @ G @ H)

G = gower(D); SST = np.trace(G)
FULL = ["sub", "trimming", "inference"] + (["aligner"] if len(ALN) > 1 else [])
ss_full = ss(G, df, FULL); ss_o = ss(G, df, ["orthology"]); ss_s = ss(G, df, ["sub"])
rng = np.random.default_rng(1); P = 9999
units = df["sub"].unique(); unit_orth = df.groupby("sub").orthology.first()
a, r_ = df.orthology.nunique(), len(REPS)
def F_orth(lbl):  # nested pseudo-F: orthology tested against subsample-within-orthology
    d = df.assign(orthology=df["sub"].map(lbl)); so = ss(G, d, ["orthology"])
    return (so / (a - 1)) / ((ss_s - so) / (a * (r_ - 1)))
obs = F_orth(unit_orth) if r_ > 1 else np.nan; cnt = 0
for _ in range(P if r_ > 1 else 0):
    cnt += F_orth(pd.Series(rng.permutation(unit_orth.values), index=unit_orth.index)) >= obs
p_orth = (cnt + 1) / (P + 1)

def perm_within(col, strata, fs_red, P=999):
    """partial R2 of `col` beyond fs_red, permuting `col` labels within strata"""
    def part(d): return (ss(G, d, fs_red + [col]) - ss(G, d, fs_red)) / SST
    o = part(df); c = 0; g = df.groupby(strata).indices
    for _ in range(P):
        v = df[col].values.copy()
        for ix in g.values(): v[ix] = v[rng.permutation(ix)]
        c += part(df.assign(**{col: v})) >= o
    return o, (c + 1) / (P + 1)
df["rep_in_orth"] = df["sub"]
R2_sub, p_sub = (np.nan, np.nan) if r_ == 1 else perm_within("rep_in_orth", ["orthology", "trimming", "inference"], ["orthology", "trimming", "inference"])
R2_inf, p_inf = perm_within("inference", ["sub", "trimming"], ["sub", "trimming"])
R2_trim, p_trim = perm_within("trimming", ["sub", "inference"], ["sub", "inference"])
vp = pd.DataFrame([
    dict(component="orthology", R2=ss_o / SST, test="nested pseudo-F vs subsamples (whole-subsample permutation)", F=obs, p_perm=p_orth),
    dict(component="subsample_within_orthology", R2=(ss_s - ss_o) / SST, test="partial, reps permuted within orth x trim x inference", F=np.nan, p_perm=p_sub),
    dict(component="trimming", R2=R2_trim, test="partial", F=np.nan, p_perm=p_trim),
    dict(component="inference", R2=R2_inf, test="partial", F=np.nan, p_perm=p_inf),
    dict(component="residual_incl_interactions", R2=1 - ss_full / SST, test="", F=np.nan, p_perm=np.nan)])
vp["marginal_R2_inference"] = np.nan
vp.loc[vp.component == "inference", "marginal_R2_inference"] = ss(G, df, ["inference"]) / SST
# interaction shares (how much of the residual is method x subsample, i.e. methods disagree differently per subsample)
df["sub_x_inf"] = df["sub"] + ":" + df.inference
df["orth_x_inf"] = df.orthology + ":" + df.inference
vp_int = pd.DataFrame([
    dict(component="orthology x inference", R2=(ss(G, df, FULL + ["orth_x_inf"]) - ss_full) / SST),
    dict(component="subsample x inference (beyond orth x inf)", R2=(ss(G, df, FULL + ["sub_x_inf"]) - ss(G, df, FULL + ["orth_x_inf"])) / SST)])
pd.concat([vp, vp_int]).to_csv(f"{OD}/variance_partition.tsv", sep="\t", index=False)

# --- per-replicate R2 (what a single-subsample study would have reported) and the resampling null
def r2_sub(ix, col):
    d = df.loc[ix].reset_index(drop=True); Gs = gower(D[np.ix_(ix, ix)])
    return ss(Gs, d, [col]) / np.trace(Gs)
per = []
for k in REPS:
    ix = np.where(df.replicate == k)[0]
    per.append(dict(replicate=k, n_trees=len(ix), R2_orthology=r2_sub(ix, "orthology"), R2_inference=r2_sub(ix, "inference"),
                    R2_trimming=r2_sub(ix, "trimming")))
per = pd.DataFrame(per); per.to_csv(f"{OD}/per_replicate_R2.tsv", sep="\t", index=False)
null = []
for o in sorted(df.orthology.unique()):
    for combo in itertools.combinations(REPS, 3):
        ix = np.where((df.orthology == o) & df.replicate.isin(combo))[0]
        null.append(dict(orthology_set=o, reps=",".join(map(str, combo)), R2_fake_orthology=r2_sub(ix, "replicate"),
                         R2_inference=r2_sub(ix, "inference")))
null = pd.DataFrame(null); null.to_csv(f"{OD}/null_orthology_R2.tsv", sep="\t", index=False)

# --- RF between trees from the same method+trimming, different subsamples
wb = []
for (t, inf), g in df.groupby(["trimming", "inference"]):
    for i, j in itertools.combinations(g.index, 2):
        wb.append(dict(trimming=t, inference=inf, same_set=df.orthology[i] == df.orthology[j], RF=D[i, j]))
wb = pd.DataFrame(wb)
wbs = wb.groupby(["inference", "same_set"]).RF.mean().unstack().reindex(columns=[False, True]).rename(columns={True: "RF_same_set_diff_subsample", False: "RF_diff_set"})
wbs.loc["ALL"] = [wb[~wb.same_set].RF.mean(), wb[wb.same_set].RF.mean()]
wbs.to_csv(f"{OD}/rf_within_between.tsv", sep="\t")

# --- clade frequencies by subsample
allc = {}
for i, cs in enumerate(C):
    for c in cs: allc.setdefault(c, []).append(i)
rec = []
for c, ix in allc.items():
    fr = len(ix) / n
    if fr < 0.05 or fr > 0.98: continue
    pres = np.zeros(n); pres[ix] = 1
    d = dict(clade=",".join(sorted(c)), freq=fr)
    for s_, g in df.groupby("sub"): d[s_] = pres[g.index].mean()
    by_o = df.assign(p=pres).groupby("orthology").p.mean(); by_s = df.assign(p=pres).groupby("sub").p.mean()
    d["effect_orthology"] = by_o.max() - by_o.min()
    d["sd_between_subsamples_within_set"] = df.assign(p=pres).groupby(["orthology", "replicate"]).p.mean().groupby(level=0).std().mean()
    d["effect_inference"] = (lambda v: v.max() - v.min())(df.assign(p=pres).groupby("inference").p.mean())
    rec.append(d)
cl = pd.DataFrame(rec).sort_values("freq", ascending=False); cl.to_csv(f"{OD}/clades_by_subsample.tsv", sep="\t", index=False)

# --- locus overlap between subsamples, via shared protein sequences
def seqset(dirn):
    s = set()
    for f in glob.glob(f"{dirn}/*.fa"):
        txt = open(f).read().split(">")[1:]
        s.add(frozenset(x.split("\n", 1)[1].replace("\n", "") for x in txt))
    return s
loc = {s_: seqset(f"results/loci/{s_}") for s_ in sorted(df["sub"].unique())}
def shared(A, B):  # loci sharing >= 1 identical protein sequence
    seqA = [set(x) for x in A]; allB = set().union(*B) if B else set()
    return sum(1 for x in seqA if x & allB)
ov = pd.DataFrame([[shared(loc[x], loc[y]) for y in loc] for x in loc], index=list(loc), columns=list(loc))
ov.to_csv(f"{OD}/locus_overlap.tsv", sep="\t")
df.drop(columns="file").to_csv(f"{OD}/trees.tsv", sep="\t", index=False)

# --- figures
cols = {"orthology": "#4C78A8", "subsample_within_orthology": "#72B7B2", "trimming": "#B5B5B5", "inference": "#F58518", "residual_incl_interactions": "#E0E0E0"}
fig, ax = plt.subplots(figsize=(7, 2.6)); left = 0
for r in vp.itertuples():
    w = max(r.R2, 0); ax.barh(0, w, left=left, color=cols[r.component], label=f"{r.component} ({w:.0%})"); left += w
ax.set_xlim(0, 1); ax.set_yticks([]); ax.set_xlabel("Share of topological variance (RF, dbRDA)")
ax.set_title(f"{n} species trees: {a} orthology sets x {r_} subsamples x {len(TRIM)} trimming x {df.inference.nunique()} methods", fontsize=9)
ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.45)); fig.tight_layout()
fig.savefig(f"{OD}/fig_reps_variance.png", dpi=200)

fig, ax = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
for a_, col, lab in [(ax[0], "R2_orthology", "orthology (real sets)"), (ax[1], "R2_inference", "inference method")]:
    a_.scatter(np.zeros(len(per)) + 0, per[col], color="#4C78A8", zorder=3, label="each replicate")
    if col == "R2_orthology" and len(null):
        a_.scatter(np.ones(len(null)) + rng.uniform(-.08, .08, len(null)), null.R2_fake_orthology, s=10, color="#B5B5B5", label="3 subsamples of the same set")
        a_.set_xticks([0, 1], ["real orthology", "null: same set"])
    else:
        a_.set_xticks([0], ["per replicate"])
    a_.set_title(lab, fontsize=9); a_.set_xlim(-.5, 1.5)
ax[0].set_ylabel("R² (single-subsample analysis)"); ax[0].legend(frameon=False, fontsize=7); fig.tight_layout()
fig.savefig(f"{OD}/fig_reps_R2_per_replicate.png", dpi=200)

if len(cl):
    subs = sorted(df["sub"].unique(), key=lambda s: (s.rsplit("-r", 1)[0], int(s.rsplit("-r", 1)[1])))
    M = cl[subs].values
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.32 * len(cl) + 1.6)))
    im = ax.imshow(M, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(subs)), subs, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cl)), [f"{c[:70]} ({f:.2f})" for c, f in zip(cl.clade, cl.freq)], fontsize=6)
    fig.colorbar(im, ax=ax, label="clade frequency"); fig.tight_layout(); fig.savefig(f"{OD}/fig_reps_clades.png", dpi=200)

pd.set_option("display.width", 250)
print(pd.concat([vp, vp_int]).to_string(index=False)); print(per.to_string(index=False))
if len(null): print(null.groupby("orthology_set")[["R2_fake_orthology", "R2_inference"]].describe().T.to_string())
print(wbs.to_string())
print(cl.drop(columns=[c for c in cl.columns if "-r" in c]).to_string(index=False))
