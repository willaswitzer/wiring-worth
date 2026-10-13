#!/bin/sh
# Session 8: distance-matched shuffled-wiring control (C. Schneider-Mizell's suggestion).
# Postsynaptic functional properties are permuted, per presynaptic neuron, only among candidate targets of
# the same cell type x area and the same 50 um soma-distance bin, so the null keeps the real targets'
# spatial distribution. Run on the main truncation sets (radial f = 1, 0.5, 0.3) and the calibrated
# sets (f = 1 and calibrated, reps rotated). ~1.5 h at 4 workers.
# Launch (after run_edit_history.sh has finished):
#   nohup caffeinate -i scripts/run_distance_null.sh > data/distance_null.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }
stamp "radial, distance-matched null"; TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 1,0.5,0.3
stamp "calibrated, distance-matched null"; TRUNC_NULL=distance TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib $R scripts/like_to_like_trunc.R perm poisson 200 4 1,0
stamp "done"
