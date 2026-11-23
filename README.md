# phylo-multiverse

Data and code for **"Phylogenies as Conditional Hypotheses: A Replicated Multiverse Analysis Separates Analytical Decisions from Locus Sampling in Phylogenomics"** (Aguirre-Carvajal & Armijos-Jaramillo, in preparation).

The study builds a "multiverse" of phylogenomic pipelines for 23 budding yeasts (Saccharomycotina; outgroup *Lipomyces starkeyi*). It then asks how much each analytical decision changes the species tree, and whether those effects persist once locus sampling is replicated.

| Experiment | Design | Species trees |
|---|---|---|
| 1. Full factorial | 3 orthology criteria × 3 aligners × 4 trimming strategies × 5 species-tree methods/models; one 50-locus sample per criterion | 180 |
| 2. Replicated subsamples | 3 criteria × 6 disjoint 50-locus subsamples × MAFFT × {no trimming, ClipKIT} × 5 methods | 180 |
| 3. Locus gradient | 3 criteria × {50, 100, 200, 400} loci × 5 subsamples × MAFFT × 4 methods, plus complete locus pools | 248 |

**Orthology criteria:**
- `of_strict`: OrthoFinder single-copy orthogroups present in all 23 taxa.
- `of_relaxed`: OrthoFinder orthogroups that are single-copy where present and present in at least 20 taxa.
- `busco`: BUSCO saccharomycetes_odb10 complete single-copy genes found in at least 20 taxa.

**Species-tree methods:**
- `concat_LG`, `concat_MFP`, `concat_C20`: concatenation in IQ-TREE 3 under LG+G4, the ModelFinder-selected model, or LG+C20+F+G4 (PMSF).
- `astral_LG`, `astral_MFP`: ASTRAL-IV on LG+G4 or ModelFinder gene trees.

## Repository layout

```
trees/
  exp1_full_factorial/
    species_trees/                       180 Newick trees, one per pipeline
    species_trees_concordance_gCF_sCF/   concatenated trees labelled with gCF/sCF (IQ-TREE)
    gene_trees/                          50 gene trees per file (input to ASTRAL and gCF)
    all_species_trees.nwk                all species trees, one per line
    species_trees_metadata.tsv           factor levels of each tree (same order as all_species_trees.nwk)
  exp2_replicate_subsamples/             same structure (subsample = 1..6)
  exp3_locus_gradient/                   same structure (n_loci = 50/100/200/400/all; no concordance trees)
  summary/
    exp1_majority_rule_*.nwk             majority-rule tree of Experiment 1 (codes / species names), rooted on Lipomyces
    complete_pool_concat_LG_rooted_species.nwk   tree from all strict OrthoFinder loci (identical topology)
alignments/      exp1 (raw and trimmed, all aligners) and exp2 (MAFFT, untrimmed and ClipKIT) alignments, tar.gz
data/
  taxa.tsv                    UniProt proteome IDs, 4-letter taxon codes, species, lineage
  orthology/                  OrthoFinder Orthogroups tables (run 1 = Exp. 1; run 2 = Exp. 2-3) and BUSCO full tables
  locus_sets/                 locus IDs of every sample/subsample, and pool sizes
analysis/        variance partitions, RF matrices, clade frequencies, per-tree tables for each experiment
figures/         manuscript figures 1-6
supplementary/   Table S1 (proteomes, sizes and BUSCO completeness)
workflow/        Snakemake workflows (Snakefile: Exp. 1-2; Snakefile_grad: Exp. 3), configs and scripts
run/             driver scripts used to run each experiment
environment.yml  software versions
```

**File naming:**
- Experiment 1: `{orthology}__{aligner}__{trimming}__{method}_{model}.nwk`.
- Experiment 2: `{orthology}-r{k}__…`, where `r0`-`r5` are subsamples 1-6.
- Experiment 3: `{orthology}-n{N}r{k}__{method}_{model}.nwk`, and `{orthology}-all__…` for the complete pools.
- Trees are unrooted as inferred and use the 4-letter taxon codes listed in `data/taxa.tsv`. Root on `Lsta` (*Lipomyces starkeyi*) for display.

## Reproducing the analyses

```bash
micromamba create -f environment.yml && micromamba activate phylo-multiverse
# 1. download the 23 UniProt proteomes listed in data/taxa.tsv into data/proteomes/{code}.fa
# 2. run OrthoFinder and BUSCO (see run/run_exp1_factorial.sh), then:
snakemake -s workflow/Snakefile -j 4 --resources mem_gb=14                                   # Experiment 1
snakemake -s workflow/Snakefile --configfile workflow/config_reps.yaml -j 4 --resources mem_gb=14   # Experiment 2
snakemake -s workflow/Snakefile_grad -j 4                                                    # Experiment 3
python workflow/scripts/analyze.py; python workflow/scripts/analyze_reps.py; python workflow/scripts/analyze_grad.py
```

All maximum-likelihood searches used IQ-TREE `--fast` with seed 1. The two OrthoFinder runs, on the same proteomes in two separate installations, produced different locus pools (strict 933 vs 1,027; relaxed 1,709 vs 1,932), whereas the BUSCO pool was identical. Both orthogroup tables are provided, and the exact locus IDs of every sample are in `data/locus_sets/`.

## Taxa

| Code | Species | Lineage |
|---|---|---|
| Scer | *Saccharomyces cerevisiae* | Saccharomycetaceae |
| Cgla | *Nakaseomyces glabratus* | Saccharomycetaceae |
| Klac | *Kluyveromyces lactis* | Saccharomycetaceae |
| Egos | *Eremothecium gossypii* | Saccharomycetaceae |
| Lthe | *Lachancea thermotolerans* | Saccharomycetaceae |
| Huva | *Hanseniaspora uvarum* | Saccharomycodaceae |
| Wano | *Wickerhamomyces anomalus* | Phaffomycetales |
| Kpha | *Komagataella phaffii* | Pichiaceae |
| Opar | *Ogataea parapolymorpha* | Pichiaceae |
| Kcap | *Kuraishia capsulata* | Pichiaceae |
| Calb | *Candida albicans* | CUG-Ser1 |
| Dhan | *Debaryomyces hansenii* | CUG-Ser1 |
| Ssti | *Scheffersomyces stipitis* | CUG-Ser1 |
| Clus | *Clavispora lusitaniae* | CUG-Ser1 |
| Bino | *Babjeviella inositovora* | CUG-Ser1 |
| Ptan | *Pachysolen tannophilus* | CUG-Ala |
| Arub | *Ascoidea rubescens* | CUG-Ser2 |
| Scra | *Saccharomycopsis crataegensis* | CUG-Ser2 |
| Ylip | *Yarrowia lipolytica* | Dipodascales |
| Slig | *Sugiyamaella lignohabitans* | Dipodascales |
| Tcas | *Tortispora caseinolytica* | Trigonopsidales |
| Nful | *Nadsonia fulvescens* | Dipodascales |
| Lsta | *Lipomyces starkeyi* | Lipomycetales |
