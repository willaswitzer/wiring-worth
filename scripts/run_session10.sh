#!/bin/sh
# Session 10: the runs left before the preprint (~14 h at 4 workers; every step resumes if interrupted).
#  1. Distance-matched and cell-type shuffles at the Feb 2026 history state (the state where the
#     orientation shift's CI excludes 0).
#  2. Session 6 baseline: Poisson shuffles topped up to 500 for every variant (Ld already has 500).
#  3. Co-registration match-rule sensitivity: lowest-residual multi-match rule and Ding's match rule
#     (scripts/match_sensitivity.py), real fits + 200 shuffles each (main variant, Poisson).
#  4. Ld (Tweedie) under the distance-matched shuffle: complete axons, radial f = 0.3, Nov 2025 state.
#  5. Ld (Tweedie) paired bootstrap of the truncation shift (radial 1, 0.5, 0.3 and both nulls), so
#     Figure 2's Ld bands are bootstrap CIs, not Wald. Longest step; last so it can't hold up the rest.
# Launch from the repo root:
#   nohup caffeinate -i scripts/run_session10.sh > data/session10.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }
HIST="TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist"

stamp "1a Feb 2026: distance-matched shuffle"
env $HIST TRUNC_LEVELS=202602,1822 TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 202602
stamp "1b Feb 2026: cell-type shuffle"
env $HIST TRUNC_LEVELS=202602,1822 $R scripts/like_to_like_trunc.R perm poisson 200 4 202602

for V in main soma ext_dend; do
  stamp "2 Session 6 Poisson shuffles to 500: $V"
  $R scripts/like_to_like.R perm $V poisson 500 600 4
done

stamp "3 match-rule variants"
$PY scripts/match_sensitivity.py
for T in residual ding; do
  $PY scripts/prepare_like_to_like.py $T
  stamp "3 $T: real fits";  L2L_TAG=$T $R scripts/like_to_like.R real
  stamp "3 $T: shuffles";   L2L_TAG=$T $R scripts/like_to_like.R perm main poisson 200 600 4
done

stamp "4a Ld real fits on history states"
env $HIST TRUNC_LEVELS=202511,202602,1822 $R scripts/like_to_like_trunc.R fit tweedie 3
stamp "4b Ld distance-matched shuffle: complete, f = 0.3"
TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm tweedie 200 3 1,0.3
stamp "4c Ld distance-matched shuffle: Nov 2025"
env $HIST TRUNC_LEVELS=202511,1822 TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm tweedie 200 3 202511

stamp "5 Ld paired bootstrap"
TRUNC_LEVELS=1,0.5,0.3 $R scripts/like_to_like_trunc.R boot tweedie 200 3
stamp "done"
