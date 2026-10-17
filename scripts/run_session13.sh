#!/bin/sh
# Session 13: like-to-like by presynaptic x postsynaptic cell type on complete axons (~40 min).
# Compares directly with Ding et al. 2025's V1 L2/3 -> L2/3 orientation analysis (their Extended Data Fig. 10).
# For 23P->23P, 4P->4P, 4P->23P, 23P->4P: real fits, distance-matched and cell-type shuffles (200 each).
# Launch from the repo root:
#   nohup caffeinate -i scripts/run_session13.sh > data/session13.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }

for C in 23P:23P 4P:4P 4P:23P 23P:4P; do
  PRE=${C%%:*}; POST=${C##*:}
  export TRUNC_PRE_CT=$PRE TRUNC_POST_CT=$POST TRUNC_LEVELS=1
  stamp "$PRE -> $POST: fit";                  $R scripts/like_to_like_trunc.R fit poisson 4
  stamp "$PRE -> $POST: distance-matched";     TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 1
  stamp "$PRE -> $POST: cell-type shuffle";    $R scripts/like_to_like_trunc.R perm poisson 200 4 1
done
unset TRUNC_PRE_CT TRUNC_POST_CT TRUNC_LEVELS
stamp "table"; $PY scripts/summarize_session13.py
stamp "done"
