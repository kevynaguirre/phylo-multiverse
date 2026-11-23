#!/bin/bash
# Locus-number gradient driver. Resumable (Snakemake).
set -u
cd /home/user/multiverso
export PATH=/opt/mm/root/envs/phylo/bin:$PATH
log(){ echo "[$(date -u +%H:%M:%S)] $*" >> grad.log; }
sync(){ tar czf /tmp/grad_snap.tgz grad.log workflow results/grad results/pool results/galn results/ggt results/gastral results/gctree 2>/dev/null \
        && mkdir -p /mnt/project-files/multiverso/gradiente && cp /tmp/grad_snap.tgz /mnt/project-files/multiverso/gradiente/grad_snapshot.tgz; }
( while true; do sleep 1800; sync; done ) & SYNCPID=$!
log "relaunch (phased)"
# phase 1: everything except the full-pool concatenated trees; phase 2: full-pool trees
T1=$(ls results/grad/*.txt | xargs -n1 basename | sed 's/.txt$//' | grep -v -- '-all$' | awk '{print "results/gctree/"$1"/LG.treefile results/gctree/"$1"/MFP.treefile"}')
snakemake -s workflow/Snakefile_grad -j 4 --keep-going --rerun-incomplete --rerun-triggers mtime --quiet rules -- $T1 results/gastral/of_relaxed-all/LG.tre results/gastral/of_relaxed-all/MFP.tre > snakemake_grad.log 2>&1
log "phase1 exit $?"
python workflow/scripts/analyze_grad.py > results/analysis_grad_stdout_phase1.txt 2>&1
snakemake -s workflow/Snakefile_grad -j 4 --keep-going --rerun-incomplete --rerun-triggers mtime --quiet rules > snakemake_grad2.log 2>&1
log "snakemake exit $?"
kill $SYNCPID; sync
python workflow/scripts/analyze_grad.py > results/analysis_grad_stdout.txt 2>&1
log "analysis exit $?"
mkdir -p /mnt/project-files/multiverso/gradiente && cp -r workflow run_grad.sh results/analysis_grad results/analysis_grad_stdout.txt results/grad/subsamples.tsv /mnt/project-files/multiverso/gradiente/
log "DONE"
