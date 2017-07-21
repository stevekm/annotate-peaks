#!/bin/bash

# Example script for running the annotations with qsub on HPC cluster

# parent peaks directory
peaks_results_dir="./peaks_results"

# find the sub-directories
peaks_dirs="$(find peaks_results/ -mindepth 1 -maxdepth 1 -type d ! -name ".db")"

# qsub log output location
log_dir="logs"
mkdir -p "$log_dir"

for peaks_dir in $peaks_dirs; do
    (
    echo "$peaks_dir"
    qsub -wd $PWD -o :${log_dir}/ -e :${log_dir}/ -j y -N "$(basename "$peaks_dir")" <<E0F
# set the correct R version if needed
module unload r
module load r/3.3.0
set -x

# need to search for specific .bed files because others which we dont want are also present
find "$peaks_dir" -type f -name "peaks.bed" | xargs ./biomaRt_ChIPpeakAnno/annotate.R
E0F

# sleep to avoid file locks
sleep 3
    )
done
