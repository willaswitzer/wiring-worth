#!/bin/sh
# Bring both sensitivity variants' Ld (Tweedie) models to 500 shuffles, then refresh the
# table and figure. Resumable: each run continues its perm_slopes_<variant>_tweedie.csv.
# Main run (4 workers) took 16.7 s/perm; expected ~2.2 h (soma, 480 to go) + ~1.8 h
# (ext_dend, 464 to go, smaller data). Budgets cap when new batches start (min).
# Launch:  nohup caffeinate -i scripts/run_sensitivity_ld_500.sh > data/like_to_like/overnight_sensitivity.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
$R scripts/like_to_like.R perm soma     tweedie 500 180 4
$R scripts/like_to_like.R perm ext_dend tweedie 500 150 4
.venv/bin/python scripts/summarize_like_to_like.py > /dev/null
echo "done: table and figure regenerated; 06_results.md text still needs the new sensitivity p-values"
