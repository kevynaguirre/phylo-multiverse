#!/bin/bash
# Replicate-subsample extension: Snakemake on {set}-r{k} + nested variance analysis. Resumable.
set -u
cd /home/user/multiverso
export PATH=/opt/mm/root/envs/phylo/bin:$PATH
log(){ echo "[$(date -u +%H:%M:%S)] $*" >> reps.log; }
sync(){ tar czf /tmp/reps_snap.tgz reps.log workflow results/loci results/*/*-r[0-9]* results/*/*-pilot data/of_out/*/Orthogroups/Orthogroups.tsv data/busco_out/*/run_*/full_table.tsv 2>/dev/null \
        && mkdir -p /mnt/project-files/multiverso/reps && cp /tmp/reps_snap.tgz /mnt/project-files/multiverso/reps/reps_snapshot.tgz; }
( while true; do sleep 1200; sync; done ) & SYNCPID=$!
log "snakemake start"
snakemake -s workflow/Snakefile --configfile workflow/config_reps.yaml -j 4 --resources mem_gb=14 --keep-going --rerun-incomplete --rerun-triggers mtime --quiet rules > snakemake_reps.log 2>&1
log "snakemake exit $?"
python workflow/scripts/analyze_reps.py > results/analysis_reps_stdout.txt 2>&1
log "analysis exit $?"
kill $SYNCPID; sync
mkdir -p /mnt/project-files/multiverso/reps
cp -r workflow run_reps.sh results/analysis_reps results/analysis_reps_stdout.txt results/loci/pools_reps.tsv /mnt/project-files/multiverso/reps/
log "DONE"
