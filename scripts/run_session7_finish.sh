#!/bin/sh
# Session 7, final runs on the pinned environment (~5 h at 4 workers; overnight):
#   1. Poisson bootstrap of the main truncation analysis, nulls' reps rotated (~2.7 h)
#   2. calibrated truncation (real partial-axon profiles per cell type): sets, fits (pooled, 23P, 4P),
#      paired bootstrap and shuffled control (~1 h)
#   3. radial truncation within 23P and 4P, paired bootstrap (~20 min)
#   4. real partial vs. full axons within 23P and 4P, bootstrap + shuffled control (~1 h)
#   5. figures and tables
# Every R step except 3-4 resumes if interrupted.
# Launch:  nohup caffeinate -i scripts/run_session7_finish.sh > data/session7_finish.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }

stamp "boot poisson (main)";       $R scripts/like_to_like_trunc.R boot poisson 200 4
stamp "calibrated sets";           $PY scripts/calibrated_truncation.py --reps 20 --seed 11
stamp "calib fits";                TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib $R scripts/like_to_like_trunc.R fit poisson 4
for ct in 23P 4P; do               TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib TRUNC_PRE_CT=$ct $R scripts/like_to_like_trunc.R fit poisson 4; done
stamp "calib boot";                TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib $R scripts/like_to_like_trunc.R boot poisson 200 4
stamp "calib perm";                TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib $R scripts/like_to_like_trunc.R perm poisson 100 4 1,0
stamp "radial within cell type";   $R scripts/celltype_truncation_check.R 200
stamp "partial vs full";           $PY scripts/prepare_partial.py
                                   $R scripts/partial_vs_full.R 200 100 4
stamp "figures";                   $PY scripts/summarize_truncation.py
                                   $PY scripts/figure_truncation_shift.py
stamp "done"
