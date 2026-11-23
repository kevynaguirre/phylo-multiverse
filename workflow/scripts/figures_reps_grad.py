"""Figures 5-6 and summary numbers for the replicate-subsample (Exp. 2) and locus-gradient (Exp. 3) experiments.
Usage: python figures_reps_grad.py <reps_dir_with_results> <grad_dir_with_results> <reps_analysis_dir> <grad_analysis_dir> <outdir>"""
import sys, glob, os, itertools, json
import numpy as np, pandas as pd, dendropy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                     "xtick.color": "#52514e", "ytick.color": "#52514e"})
INK, MUTED, GRID, LIGHT = "#0b0b0b", "#52514e", "#e4e3df", "#b9b8b3"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RD, GD, RA, GA, OUTD = sys.argv[1:6]
OUTG = "Lsta"
CL = {"N1a": "Arub,Cgla,Egos,Huva,Klac,Lthe,Scer,Scra,Wano", "N1b": "Arub,Bino,Calb,Clus,Dhan,Kcap,Kpha,Opar,Ptan,Scra,Ssti",
      "N2a": "Nful,Ylip", "N2b": "Nful,Slig", "N3a": "Egos,Klac,Lthe", "N3b": "Cgla,Lthe,Scer"}
OLAB = {"of_strict": "OrthoFinder strict", "of_relaxed": "OrthoFinder relaxed", "busco": "BUSCO"}

def load(files):
    tns = dendropy.TaxonNamespace(); out = []
    for f in files:
        s = open(f).read().strip()
        out.append(dendropy.Tree.get(data=s, schema="newick", taxon_namespace=tns, preserve_underscores=True))
    taxa = sorted(t.label for t in tns)
    res = []
    for tr in out:
        tr.encode_bipartitions(); cs = set()
        for b in tr.bipartition_encoding:
            side = frozenset(t.label for t in b.leafset_taxa(tns))
            if OUTG in side: side = frozenset(taxa) - side
            if 1 < len(side) < len(taxa) - 1: cs.add(",".join(sorted(side)))
        res.append(cs)
    return res

# ---------------- Experiment 2: replicate subsamples ----------------
rows = []
for f in glob.glob(f"{RD}/results/ctree/*-r[0-5]/mafft/*/*.treefile"):
    s, a, t, m = f.split("/")[-4:]; rows.append(dict(sub=s, trimming=t, inference="concat_" + m[:-9], file=f))
for f in glob.glob(f"{RD}/results/astral/*-r[0-5]/mafft/*/*.tre"):
    s, a, t, m = f.split("/")[-4:]; rows.append(dict(sub=s, trimming=t, inference="astral_" + m[:-4], file=f))
R = pd.DataFrame(rows); R = R[R.trimming.isin(["none", "clipkit"])].reset_index(drop=True)
R["orthology"] = R["sub"].str.rsplit("-r", n=1).str[0]; R["rep"] = R["sub"].str.rsplit("-r", n=1).str[1].astype(int)
R["method"] = R.inference.str.split("_").str[0]
C = load(R.file)
for k, c in CL.items(): R[k] = [c in cs for cs in C]
assert len(R) == 180, len(R)
S = {}
S["reps_ntrees"] = len(R)
S["reps_N1a_by_sub"] = R.groupby("sub").N1a.mean().round(2).to_dict()
S["reps_N1a_by_sub_method"] = R.groupby(["sub", "method"]).N1a.mean().unstack().round(2).to_dict()
S["reps_N1a_overall"] = R.N1a.mean(); S["reps_N1a_by_orth"] = R.groupby("orthology").N1a.mean().round(3).to_dict()
S["reps_N2a_by_method"] = R.groupby("method").N2a.mean().round(3).to_dict()
S["reps_N2a_by_orth_method"] = R.groupby(["orthology", "method"]).N2a.mean().round(3).unstack().to_dict()
S["reps_N2a_by_sub_method"] = R.groupby(["sub", "method"]).N2a.mean().unstack().round(2).to_dict()
S["reps_N3a"] = R.N3a.mean()
# how many subsamples per set give majority N1a / N1b
bysub = R.groupby(["orthology", "sub"]).N1a.mean()
S["reps_subs_N1a_majority"] = {o: int((g > 0.5).sum()) for o, g in bysub.groupby(level=0)}
S["reps_subs_N1_mixed"] = {o: int(((g > 0) & (g < 1)).sum()) for o, g in bysub.groupby(level=0)}
vp = pd.read_csv(f"{RA}/variance_partition.tsv", sep="\t").set_index("component")
per = pd.read_csv(f"{RA}/per_replicate_R2.tsv", sep="\t")
null = pd.read_csv(f"{RA}/null_orthology_R2.tsv", sep="\t")
S["per_rep_R2_orth"] = per.R2_orthology.round(3).tolist(); S["null_R2_mean"] = null.R2_fake_orthology.mean()
S["null_R2_range"] = [null.R2_fake_orthology.min(), null.R2_fake_orthology.max()]
S["null_R2_by_set"] = null.groupby("orthology_set").R2_fake_orthology.mean().round(3).to_dict()
S["null_frac_ge_0499"] = float((null.R2_fake_orthology >= 0.499).mean())

# ---------------- Experiment 3: locus gradient ----------------
rows = []
for f in glob.glob(f"{GD}/results/gctree/*/*.treefile"):
    s = f.split("/")[-2]; rows.append(dict(sub=s, inference="concat_" + os.path.basename(f)[:-9], file=f))
for f in glob.glob(f"{GD}/results/gastral/*/*.tre"):
    s = f.split("/")[-2]; rows.append(dict(sub=s, inference="astral_" + os.path.basename(f)[:-4], file=f))
Gt = pd.DataFrame(rows)
Gt = Gt[Gt.file.map(lambda f: os.path.getsize(f) > 0)].reset_index(drop=True)
meta = pd.read_csv(f"{GA}/trees.tsv", sep="\t")[["sub", "orthology", "N", "replicate"]].drop_duplicates("sub").set_index("sub")
Gt = Gt[Gt["sub"].isin(meta.index)].reset_index(drop=True)
Gt["orthology"] = Gt["sub"].map(meta.orthology); Gt["N"] = Gt["sub"].map(meta.N).astype(str)
Gt["method"] = Gt.inference.str.split("_").str[0]
C = load(Gt.file)
for k, c in CL.items(): Gt[k] = [c in cs for cs in C]
S["grad_ntrees"] = len(Gt); S["grad_full_trees"] = Gt[Gt.N == "all"][["orthology", "inference"]].values.tolist()
S["grad_N1a"] = Gt.groupby(["N", "method"]).N1a.mean().round(3).unstack().to_dict()
S["grad_N2a"] = Gt.groupby(["N", "method"]).N2a.mean().round(3).unstack().to_dict()
S["grad_N3a"] = Gt.groupby(["N", "method"]).N3a.mean().round(3).unstack().to_dict()
S["grad_full_clades"] = Gt[Gt.N == "all"][["N1a", "N2a", "N3a"]].all().to_dict()
vN = pd.read_csv(f"{GA}/variance_by_N.tsv", sep="\t"); rfN = pd.read_csv(f"{GA}/rf_by_N.tsv", sep="\t")
S["vN"] = vN.round(3).to_dict(orient="records")
json.dump(S, open(f"{OUTD}/reps_grad_summary.json", "w"), indent=1, default=float)

# ---------------- Figure 5 ----------------
fig = plt.figure(figsize=(7.0, 5.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], width_ratios=[1.35, 1], hspace=0.62, wspace=0.42)
a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, :])
comp = [("Orthology criterion", 49.9, 100 * vp.loc["orthology", "R2"]),
        ("Locus subsample\n(within criterion)", None, 100 * vp.loc["subsample_within_orthology", "R2"]),
        ("Species-tree method\n& model", 30.4, 100 * max(vp.loc["inference", "R2"], 0)),
        ("Subsample × method", None, 100 * vp.loc["subsample x inference (beyond orth x inf)", "R2"]),
        ("Trimming", 0.4, 100 * max(vp.loc["trimming", "R2"], 0))]
y = np.arange(len(comp))[::-1]; h = 0.36
for yi, (lab, e1, e2) in zip(y, comp):
    if e1 is not None:
        a.barh(yi + h / 2 + 0.02, e1, height=h, color=LIGHT); a.text(e1 + 1.2, yi + h / 2 + 0.02, f"{e1:.1f}%", va="center", fontsize=6.6, color=MUTED)
    else:
        a.text(1.2, yi + h / 2 + 0.02, "not estimable", va="center", fontsize=6.6, color=MUTED, style="italic")
    a.barh(yi - h / 2 - 0.02, e2, height=h, color=BLUE); a.text(max(e2, 0) + 1.2, yi - h / 2 - 0.02, f"{e2:.1f}%", va="center", fontsize=6.6, color=INK)
a.set_yticks(y, [x[0] for x in comp]); a.set_xlim(0, 75); a.set_xlabel("Topological variance explained (%)")
a.xaxis.grid(True, color=GRID, lw=0.5); a.set_axisbelow(True)
a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=LIGHT), plt.Rectangle((0, 0), 1, 1, color=BLUE)],
         labels=["Exp. 1: one subsample per criterion", "Exp. 2: six subsamples per criterion"],
         frameon=False, fontsize=6.6, loc="lower right", handlelength=1)
a.set_title("A", loc="left", weight="bold", fontsize=10)
rng = np.random.default_rng(3)
b.scatter(np.zeros(len(per)) + rng.uniform(-.06, .06, len(per)), 100 * per.R2_orthology, s=22, color=BLUE, zorder=3, edgecolor="white", lw=0.5)
b.scatter(np.ones(len(null)) + rng.uniform(-.14, .14, len(null)), 100 * null.R2_fake_orthology, s=9, color=LIGHT, zorder=2)
b.axhline(49.9, color=MUTED, lw=0.7, ls=(0, (3, 2))); b.text(-0.47, 47.5, "Exp. 1\nestimate", fontsize=6.2, color=MUTED, ha="left", va="top")
for xi, v in [(0, per.R2_orthology), (1, null.R2_fake_orthology)]:
    b.plot([xi - .25, xi + .25], [100 * v.mean()] * 2, color=INK, lw=1.2, zorder=4)
b.set_xticks([0, 1], ["Three criteria,\none subsample each\n(n = 6)", "Three subsamples\nof one criterion\n(null, n = 60)"])
b.set_xlim(-.5, 1.5); b.set_ylim(0, 105); b.set_ylabel("Variance attributed to\n'orthology' (%)")
b.yaxis.grid(True, color=GRID, lw=0.5); b.set_axisbelow(True)
b.set_title("B", loc="left", weight="bold", fontsize=10)
order = ["of_strict", "of_relaxed", "busco"]; xs, lab, xt = [], [], []
x = 0
for o in order:
    for r in range(6):
        s = f"{o}-r{r}"; g = R[R["sub"] == s]
        fc = g[g.method == "concat"].N1a.mean() * 100; fa = g[g.method == "astral"].N1a.mean() * 100
        c.bar(x - 0.2, fc, width=0.38, color=BLUE); c.bar(x + 0.2, fa, width=0.38, color=ORANGE)
        for dx, v in [(-0.2, fc), (0.2, fa)]:
            if v == 0: c.text(x + dx, 2, "0", ha="center", va="bottom", fontsize=6, color=MUTED)
        xt.append(x); lab.append(f"{r + 1}"); x += 1
    c.text(x - 3.5, 112, OLAB[o], ha="center", fontsize=7.2, color=INK)
    x += 0.8
c.set_xticks(xt, lab, fontsize=6.6); c.set_xlabel("Independent 50-locus subsample")
c.set_ylim(0, 118); c.set_yticks([0, 25, 50, 75, 100]); c.set_ylabel("Pipelines recovering N1a (%)")
c.yaxis.grid(True, color=GRID, lw=0.5); c.set_axisbelow(True)
c.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)],
         labels=["Concatenation (6 pipelines)", "ASTRAL (4 pipelines)"], frameon=False, fontsize=6.6,
         loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=2, handlelength=1)
c.set_title("C", loc="left", weight="bold", fontsize=10, x=-0.06)
fig.savefig(f"{OUTD}/Fig5_replicates.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# ---------------- Figure 6 ----------------
fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.5), gridspec_kw={"wspace": 0.5})
Ns = [50, 100, 200, 400]; xp = np.arange(4)
for meth, col, lab in [("concat", BLUE, "Concatenation"), ("astral", ORANGE, "ASTRAL")]:
    v = rfN[rfN.inference.str.startswith(meth)].groupby("N").RF_between_replicates_same_set.mean().loc[Ns]
    a.plot(xp, v.values, color=col, lw=2, marker="o", ms=4.5, mec="white", mew=0.6)
    a.plot([3, 4], [v.values[-1], 0], color=col, lw=1, ls=(0, (2, 2)))
    if meth == "astral": a.text(1.1, 0.052, lab, color=INK, fontsize=6.6, va="center")
    else: a.text(1.1, 0.0125, lab, color=INK, fontsize=6.6, va="center")
a.scatter([4], [0], s=26, color=INK, zorder=4)
a.set_xticks(list(xp) + [4], [str(n) for n in Ns] + ["All"]); a.set_xlabel("Loci per subsample")
a.set_ylabel("Mean nRF between\nsubsamples of one criterion"); a.set_ylim(-0.003, 0.09)
a.yaxis.grid(True, color=GRID, lw=0.5); a.set_axisbelow(True); a.set_title("A", loc="left", weight="bold", fontsize=10)
vv = vN.set_index("N").loc[Ns]
for col_, colr, lab, mk in [("R2_subsample_within_orth", BLUE, "Locus subsample", "o"),
                             ("R2_orthology", ORANGE, "Orthology", "s"),
                             ("R2_inference_partial", AQUA, "Method", "^")]:
    yv = 100 * vv[col_].clip(lower=0).values
    b.plot(xp, yv, color=colr, lw=2, marker=mk, ms=4.5, mec="white", mew=0.6, label=lab)
b.set_xticks(xp, [str(n) for n in Ns]); b.set_xlabel("Loci per subsample"); b.set_ylabel("Topological variance\nexplained (%)")
b.set_ylim(0, 70); b.yaxis.grid(True, color=GRID, lw=0.5); b.set_axisbelow(True)
b.legend(frameon=False, fontsize=6.4, loc="upper right", handlelength=1.6)
b.set_title("B", loc="left", weight="bold", fontsize=10)
for key, colr, lab, mk in [("N1a", BLUE, "N1a (CUG-Ser2)", "o"), ("N2a", ORANGE, "N2a (Nadsonia)", "s")]:
    v = Gt.groupby("N")[key].mean()
    yv = [100 * v.loc[str(n)] for n in Ns] + [100 * v.loc["all"]]
    c.plot(list(xp) + [4], yv, color=colr, lw=2, marker=mk, ms=4.5, mec="white", mew=0.6, label=lab)
c.set_xticks(list(xp) + [4], [str(n) for n in Ns] + ["All"]); c.set_xlabel("Loci per subsample")
c.set_ylabel("Trees recovering\nthe resolution (%)"); c.set_ylim(50, 102)
c.yaxis.grid(True, color=GRID, lw=0.5); c.set_axisbelow(True); c.legend(frameon=False, fontsize=6.4, loc="lower right", handlelength=1.6)
c.set_title("C", loc="left", weight="bold", fontsize=10)
fig.savefig(f"{OUTD}/Fig6_gradient.png", dpi=600, bbox_inches="tight"); plt.close(fig)
print(json.dumps(S, indent=1, default=float))
