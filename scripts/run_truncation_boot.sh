#!/bin/sh
# Session 7: rerun the Poisson bootstrap with the nulls' reps rotated across bootstrap samples
# (sample b uses rep (b-1) mod 5), then refresh both figures. ~2.7 h at 4 workers.
# The previous bootstrap (rep 0 of the nulls only) is kept in data/truncation/run2_boot_rep0/.
# Launch:  nohup caffeinate -i scripts/run_truncation_boot.sh > data/truncation_boot.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
echo "=== $(date '+%H:%M:%S') boot poisson"
.r-env/bin/Rscript scripts/like_to_like_trunc.R boot poisson 200 4
.venv/bin/python scripts/summarize_truncation.py
.venv/bin/python scripts/figure_truncation_shift.py
echo "=== $(date '+%H:%M:%S') done"
