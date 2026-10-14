#!/bin/sh
# Session 8: uncertainty for the edit-history result (~3 h at 4 workers).
# Time points: Apr 2023 (start of the partial plateau), Nov 2025 (end of the plateau), Feb 2026 (mid-
# completion), v1822 (complete). For each:
#   1. paired bootstrap over presynaptic neurons of the shift vs the complete state (200 samples)
#   2. distance-matched shuffled-wiring control (200 per time point)
#   3. cell-type x area shuffled-wiring control (200 per time point; the Session 6 null)
# Launch:  nohup caffeinate -i scripts/run_history_ci.sh > data/history_ci.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
LV=202304,202511,202602,1822
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }
stamp "boot";              TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist TRUNC_LEVELS=$LV $R scripts/like_to_like_trunc.R boot poisson 200 4
stamp "distance null";     TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist TRUNC_LEVELS=$LV TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 202304,202511,1822
stamp "cell-type null";    TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist TRUNC_LEVELS=$LV $R scripts/like_to_like_trunc.R perm poisson 200 4 202304,202511
stamp "done"
