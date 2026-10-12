#!/bin/sh
# Session 7 sensitivity runs (~3 h at 4 workers; overnight):
#   1. calibration uncertainty: 20 calibrated-truncation reps, each with its own target profile drawn by
#      resampling real partial and full 23P/4P axons; fits + paired bootstrap (reps rotated) -> the CI on
#      the realistic bias including uncertainty in the calibration itself (~45 min)
#   2. block-size sensitivity of the random-cable-loss null: 10 um and 50 um cubes (main run: 25 um);
#      fits + paired bootstrap (~1 h each)
# Every R step resumes if interrupted.
# Launch:  nohup caffeinate -i scripts/run_session7_sensitivity.sh > data/session7_sensitivity.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }

stamp "calibration-uncertainty sets"; $PY scripts/calibrated_truncation.py --reps 20 --seed 12 --profile-boot
for t in calibunc blocks10 blocks50; do
  if [ "$t" = blocks10 ]; then stamp "blocks10 sets"; $PY scripts/truncate_like_to_like.py --reps 5 --block-um 10 --seed 8 --only-blocks blocks10; fi
  if [ "$t" = blocks50 ]; then stamp "blocks50 sets"; $PY scripts/truncate_like_to_like.py --reps 5 --block-um 50 --seed 9 --only-blocks blocks50; fi
  stamp "$t fit";  TRUNC_INDEX=index_$t.csv TRUNC_TAG=$t $R scripts/like_to_like_trunc.R fit poisson 4
  stamp "$t boot"; TRUNC_INDEX=index_$t.csv TRUNC_TAG=$t $R scripts/like_to_like_trunc.R boot poisson 200 4
done
stamp "done"
