#!/bin/sh
# Session 7 overnight run: the whole truncation pipeline, most important steps first.
# Every R step resumes from its CSV, so if it is still running when you need the Mac,
# Ctrl-C / kill it and rerun the same command later: finished work is kept.
# Expected (8 GB Mac, nothing else heavy running):
#   sets + all fits ~15 min; Poisson shuffle 200 x 3 levels ~45 min;
#   Poisson bootstrap 200 ~3.5 h; Ld (Tweedie) shuffle 200 x 3 levels ~3 h. Total ~7-7.5 h.
# The table and figure are regenerated after each stage, so partial runs are usable.
# Launch:  nohup caffeinate -i scripts/run_truncation_overnight.sh > data/truncation_overnight.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }

stamp "truncation sets";      $PY scripts/truncate_like_to_like.py --reps 5 --block-um 25 --seed 7
stamp "fit poisson";          $R scripts/like_to_like_trunc.R fit poisson 4
stamp "fit tweedie";          $R scripts/like_to_like_trunc.R fit tweedie 4
$PY scripts/summarize_truncation.py
stamp "perm poisson";         $R scripts/like_to_like_trunc.R perm poisson 200 4 1,0.5,0.3
$PY scripts/summarize_truncation.py
stamp "boot poisson";         $R scripts/like_to_like_trunc.R boot poisson 200 4
$PY scripts/summarize_truncation.py
stamp "perm tweedie";         $R scripts/like_to_like_trunc.R perm tweedie 200 4 1,0.5,0.3
$PY scripts/summarize_truncation.py
stamp "done"
