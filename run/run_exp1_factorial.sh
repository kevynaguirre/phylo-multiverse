#!/bin/bash
# Driver for the cloud pilot: OrthoFinder -> BUSCO -> locus selection -> Snakemake multiverse -> analysis
set -u
cd /home/user/multiverso
export PATH=/opt/mm/root/envs/phylo/bin:$PATH
log(){ echo "[$(date -u +%H:%M:%S)] $*" >> pilot.log; }
log "waiting for OrthoFinder"
while pgrep -f "orthofinder -f proteomes" >/dev/null; do sleep 30; done
log "OrthoFinder done: $(ls data/of_out/*/Orthogroups/Orthogroups.tsv)"
log "BUSCO start"
cd data
ls proteomes/*.fa | xargs -P 4 -I{} sh -c 'c=$(basename {} .fa); [ -s busco_out/$c/run_saccharomycetes_odb10/full_table.tsv ] || busco -i {} -m prot -l busco_db/saccharomycetes_odb10 --offline -c 1 -o $c --out_path busco_out -q > /dev/null 2>&1'
cd ..
log "BUSCO done: $(ls data/busco_out/*/run_*/full_table.tsv | wc -l) tables"
[ -s results/loci/pools.tsv ] || python workflow/scripts/select_loci.py data/proteomes data/of_out data/busco_out >> pilot.log 2>&1
log "snakemake start"
snakemake -s workflow/Snakefile -j 4 --resources mem_gb=14 --keep-going --rerun-incomplete --rerun-triggers mtime --quiet rules > snakemake.log 2>&1
log "snakemake exit $?"
python workflow/scripts/analyze.py > results/analysis_stdout.txt 2>&1
log "analysis exit $?"
mkdir -p /mnt/project-files/multiverso
cp -r workflow run_pilot.sh data/taxa.tsv /mnt/project-files/multiverso/
cp -r results/analysis results/loci/pools.tsv results/analysis_stdout.txt /mnt/project-files/multiverso/ 2>/dev/null
mkdir -p /mnt/project-files/multiverso/species_trees && tar czf /mnt/project-files/multiverso/species_trees/species_trees.tar.gz results/ctree results/astral
log "copied to project files; DONE"
