"""Publication figures for the multiverse manuscript (static, print)."""
import numpy as np, pandas as pd, dendropy, glob, yaml
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                     "xtick.color": "#52514e", "ytick.color": "#52514e"})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
OUT = "results/figures"
cfg = yaml.safe_load(open("workflow/config.yaml"))
vp = pd.read_csv("results/analysis/variance_partition.tsv", sep="\t")
tr = pd.read_csv("results/analysis/trees.tsv", sep="\t")
cl = pd.read_csv("results/analysis/clades.tsv", sep="\t")
OLAB = {"of_strict": "OrthoFinder\nstrict", "of_relaxed": "OrthoFinder\nrelaxed", "busco": "BUSCO"}
ILAB = {"concat_LG": "Concat.\nLG+G4", "concat_MFP": "Concat.\nModelFinder", "concat_C20": "Concat.\nC20 (PMSF)",
        "astral_LG": "ASTRAL\nLG+G4", "astral_MFP": "ASTRAL\nModelFinder"}
FLAB = {"orthology": "Orthology", "aligner": "Aligner", "trimming": "Trimming", "inference": "Species-tree\nmethod & model"}

# ---------- Figure 1: study design ----------
fig, ax = plt.subplots(figsize=(7.0, 2.6)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 4)
stages = [("23 reference\nproteomes", ["UniProt", "Saccharomycotina"]),
          ("Orthology", ["OrthoFinder strict", "OrthoFinder relaxed", "BUSCO"]),
          ("Alignment", ["MAFFT L-INS-i", "MUSCLE5", "FAMSA"]),
          ("Trimming", ["None", "trimAl gappyout", "ClipKIT smart-gap", "BMGE"]),
          ("Species tree\n& model", ["Concat. LG+G4", "Concat. ModelFinder", "Concat. C20 (PMSF)", "ASTRAL LG+G4", "ASTRAL ModelFinder"])]
W = [1.55, 1.75, 1.6, 1.85, 2.05]; gap = (10 - sum(W)) / 4
for i, (title, items) in enumerate(stages):
    w = W[i]; x = sum(W[:i]) + i * gap
    ax.add_patch(FancyBboxPatch((x, 0.15), w, 3.3, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc="#f4f3f0" if i else "#ffffff", ec=MUTED, lw=0.6))
    ax.text(x + w / 2, 3.15, title, ha="center", va="center", fontsize=8.5, weight="bold", color=INK)
    for j, it in enumerate(items):
        ax.text(x + w / 2, 2.5 - j * 0.42, it, ha="center", va="center", fontsize=6.9, color=INK)
    if i: ax.text(x + w / 2, 0.35, f"{len(items)} levels", ha="center", fontsize=6.8, color=MUTED, style="italic")
    if i < 4: ax.annotate("", xy=(x + w + gap - 0.03, 1.8), xytext=(x + w + 0.03, 1.8),
                          arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.8))
ax.text(10, -0.15, "3 × 3 × 4 × 5 = 180 species trees (50 loci per orthology set; 3,600 gene trees)",
        ha="right", va="top", fontsize=7.2, color=MUTED)
fig.savefig(f"{OUT}/Fig1_design.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# ---------- Figure 2: variance partition + recovery of the majority topology ----------
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1, 1.25]})
v = vp[vp.factor != "all_factors"].copy()
v["pct"] = 100 * v.partial_R2.clip(lower=0)
order = ["orthology", "inference", "trimming", "aligner"]; v = v.set_index("factor").loc[order]
y = np.arange(len(v))[::-1]
a.barh(y, v.pct, height=0.55, color=BLUE)
for yi, (f, r) in zip(y, v.iterrows()):
    p = "p = 0.001" if r.p_perm <= 0.001 else f"p = {r.p_perm:.2f}"
    a.text(r.pct + 1.2, yi, f"{r.pct:.1f}%  ({p})", va="center", fontsize=7, color=INK)
a.set_yticks(y, [FLAB[f] for f in order]); a.set_xlim(0, 75)
a.set_xlabel("Topological variance explained (%)"); a.xaxis.grid(True, color=GRID, lw=0.5); a.set_axisbelow(True)
a.set_title("A", loc="left", weight="bold", fontsize=10)
M = (tr.assign(ok=(tr.RF_to_majority == 0)).groupby(["orthology", "inference"]).ok.mean().unstack() * 100)
M = M.loc[["of_strict", "of_relaxed", "busco"], ["concat_LG", "concat_MFP", "concat_C20", "astral_LG", "astral_MFP"]]
im = b.imshow(M.values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        val = M.values[i, j]; b.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=7.5,
                                     color="white" if val > 55 else INK)
b.set_xticks(range(5), ["Concat.\nLG+G4", "Concat.\nMF", "Concat.\nC20", "ASTRAL\nLG+G4", "ASTRAL\nMF"], fontsize=6.8); b.set_yticks(range(3), [OLAB[o].replace("\n", " ") for o in M.index])
b.tick_params(length=0); [s.set_visible(False) for s in b.spines.values()]
b.set_title("B  Pipelines recovering the majority topology", loc="left", weight="bold", fontsize=9)
fig.tight_layout(); fig.savefig(f"{OUT}/Fig2_variance.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# ---------- Figure 3: decision-sensitive clades across the multiverse ----------
CN = {"Arub,Cgla,Egos,Huva,Klac,Lthe,Scer,Scra,Wano": "N1a  CUG-Ser2 + Saccharomycetaceae/Hanseniaspora/Wickerhamomyces",
      "Arub,Bino,Calb,Clus,Dhan,Kcap,Kpha,Opar,Ptan,Scra,Ssti": "N1b  CUG-Ser2 + CUG-Ser1/Pichiaceae/CUG-Ala",
      "Nful,Ylip": "N2a  Nadsonia + Yarrowia", "Nful,Slig": "N2b  Nadsonia + Sugiyamaella",
      "Egos,Klac,Lthe": "N3a  Eremothecium + Kluyveromyces + Lachancea",
      "Cgla,Lthe,Scer": "N3b  Lachancea + Saccharomyces + Nakaseomyces"}
tns = dendropy.TaxonNamespace(); rows = []
OUTG = cfg["outgroup"]
for r in tr.itertuples():
    f = (f"results/ctree/{r.orthology}/{r.aligner}/{r.trimming}/{r.model}.treefile" if r.method == "concat"
         else f"results/astral/{r.orthology}/{r.aligner}/{r.trimming}/{r.model}.tre")
    t = dendropy.Tree.get(path=f, schema="newick", taxon_namespace=tns, preserve_underscores=True); t.encode_bipartitions()
    taxa = frozenset(x.label for x in tns); cs = set()
    for bp in t.bipartition_encoding:
        side = frozenset(x.label for x in bp.leafset_taxa(tns))
        if OUTG in side: side = taxa - side
        cs.add(",".join(sorted(side)))
    rows.append([1 if k in cs else 0 for k in CN])
P = np.array(rows).T
iorder = ["concat_LG", "concat_MFP", "concat_C20", "astral_LG", "astral_MFP"]
tr["_o"] = tr.orthology.map({"of_strict": 0, "of_relaxed": 1, "busco": 2}); tr["_i"] = tr.inference.map({k: i for i, k in enumerate(iorder)})
idx = tr.sort_values(["_o", "_i", "aligner", "trimming"]).index.values; P = P[:, idx]
fig = plt.figure(figsize=(7.0, 3.1)); ax = fig.add_axes([0.43, 0.2, 0.55, 0.66])
cmap = matplotlib.colors.ListedColormap(["#f1f0ec", "#1f5fa8"])
ax.imshow(P, aspect="auto", cmap=cmap, interpolation="nearest", vmin=0, vmax=1)
for yy in [1.5, 3.5]: ax.axhline(yy, color="white", lw=3)
for xx in [59.5, 119.5]: ax.axvline(xx, color="white", lw=3)
for xx in np.arange(11.5, 180, 12): ax.axvline(xx, color="white", lw=0.8)
freqs = P.mean(1)
ax.set_yticks(range(6), [f"{lab}  ({fr*100:.0f}%)".replace("Saccharomycetaceae/", "Saccharomycetaceae/\n") for lab, fr in zip(CN.values(), freqs)], fontsize=6.4)
ax.set_xticks([]); ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
for k, o in enumerate(["of_strict", "of_relaxed", "busco"]):
    ax.text(30 + 60 * k, -1.55, OLAB[o].replace("\n", " "), ha="center", fontsize=8, weight="bold", color=INK)
    for j, inf in enumerate(iorder):
        ax.text(60 * k + 12 * j + 6, -0.75, ILAB[inf].replace("Concat.", "C").replace("ASTRAL", "A").replace("ModelFinder", "MF").replace("C20 (PMSF)", "C20").replace("LG+G4", "LG"),
                ha="center", fontsize=5.6, color=MUTED, linespacing=0.9)
ax.set_xlabel("180 pipelines, grouped by orthology and species-tree method\n(C, concatenation; A, ASTRAL; MF, ModelFinder); within each block:\naligner (MAFFT, MUSCLE5, FAMSA) × trimming (none, gappyout, ClipKIT, BMGE)", fontsize=6.4, color=MUTED)
fig.savefig(f"{OUT}/Fig3_sensitive_nodes.png", dpi=600, bbox_inches="tight"); plt.close(fig)

# ---------- Figure 4: support vs stability ----------
SENS = {k.split("  ")[0]: v for v, k in zip(CN.keys(), CN.values())}
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8))
from scipy.stats import spearmanr
for a_, col, lab in zip(axs, ["mean_gCF", "mean_localPP"], ["Mean gene concordance factor, gCF (%)", "Mean ASTRAL local posterior probability"]):
    s = cl.dropna(subset=[col]); sens = s.freq < 0.98
    a_.scatter(s[col][~sens], 100 * s.freq[~sens], s=22, color=BLUE, edgecolor="white", lw=0.8, label="Stable clade (≥ 98%)", zorder=3)
    a_.scatter(s[col][sens], 100 * s.freq[sens], s=26, color=ORANGE, edgecolor="white", lw=0.8, marker="D", label="Decision-sensitive clade", zorder=3)
    for r in s[sens].itertuples():
        nm = [k.split("  ")[0] for v, k in CN.items() if v == r.clade]
        off = {("mean_localPP", "N1a"): (-20, 5), ("mean_localPP", "N1b"): (-22, -9), ("mean_localPP", "N3b"): (6, 2)}.get((col, nm[0] if nm else ""), (5, 3))
        if nm: a_.annotate(nm[0], (getattr(r, col), 100 * r.freq), xytext=off, textcoords="offset points", fontsize=6.5, color=INK)
    rho, p = spearmanr(s.freq, s[col])
    a_.text(1.0, 1.03, f"Spearman ρ = {rho:.2f}, p = {p:.1g}, n = {len(s)} clades", transform=a_.transAxes, ha="right", fontsize=6.8, color=MUTED)
    a_.set_xlabel(lab); a_.set_ylabel("Frequency across multiverse (%)"); a_.yaxis.grid(True, color=GRID, lw=0.5); a_.set_axisbelow(True)
axs[0].legend(frameon=False, fontsize=6.8, loc="center right"); [a_.set_ylim(0, 105) for a_ in axs]
axs[0].set_title("A", loc="left", weight="bold", fontsize=10); axs[1].set_title("B", loc="left", weight="bold", fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/Fig4_support_vs_stability.png", dpi=600, bbox_inches="tight"); plt.close(fig)
print("ok")
