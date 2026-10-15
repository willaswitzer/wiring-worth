#!/bin/sh
# Session 11: follow-ups to C. Schneider-Mizell's reply (~4 h after Session 10 finishes).
#  1. Newest public release vs v1822: status of the 193 axons, newly complete 4P axons (network only; runs now).
#  2. Truncation at one absolute radius for every axon: same path radius from the soma, and a cylinder
#     around the proofread column (scripts/absolute_truncation.py; light, runs now).
#  3. Waits for scripts/run_session10.sh to finish, then: fits, cell-type and distance-matched shuffles
#     (200 each at two levels per condition), paired bootstrap (200), summary table.
#  4. Regenerates Figures 1 and 2 (Figure 1 picks up the Feb 2026 state from Session 10).
# Launch from the repo root (safe to start while Session 10 is still running):
#   nohup caffeinate -i scripts/run_session11.sh > data/session11.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }

stamp "1 newest release check"
$PY scripts/check_new_release.py || echo "release check failed (network/token?); rerun it on its own later"

stamp "2 absolute-radius truncation sets"
$PY scripts/absolute_truncation.py

stamp "waiting for run_session10.sh to finish"
while pgrep -f run_session10.sh > /dev/null; do sleep 300; done

for C in pathabs cyl; do
  if [ $C = pathabs ]; then LV=250,400; else LV=100,200; fi
  stamp "3 $C: fits";                 TRUNC_INDEX=index_$C.csv TRUNC_TAG=$C $R scripts/like_to_like_trunc.R fit poisson 4
  stamp "3 $C: distance-matched";     TRUNC_INDEX=index_$C.csv TRUNC_TAG=$C TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 $LV
  stamp "3 $C: cell-type shuffle";    TRUNC_INDEX=index_$C.csv TRUNC_TAG=$C $R scripts/like_to_like_trunc.R perm poisson 200 4 $LV
  stamp "3 $C: bootstrap";            TRUNC_INDEX=index_$C.csv TRUNC_TAG=$C $R scripts/like_to_like_trunc.R boot poisson 200 4
done
$PY scripts/summarize_abs_truncation.py

stamp "4 figures"
$PY scripts/figure_history.py
$PY scripts/figure_truncation_shift.py
stamp "done"
