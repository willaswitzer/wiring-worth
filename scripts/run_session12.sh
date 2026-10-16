#!/bin/sh
# Session 12: why the orientation rule looked weaker in the partial data, and the postsynaptic rollback (~3 h).
#  1. Missing vs. attached synapses at each history state (scripts/missing_synapses.py; ~2 min).
#  2. Postsynaptic rollback: fetch postsynaptic supervoxel roots at 4 states (CAVE, ~5 min), build histpp sets.
#  3. Count-matched random loss: each complete axon loses synapses at random down to its history count (20 reps).
#  4. Fits: histsub, histpp, and the history per presynaptic cell type (4P, 23P).
#  5. Shuffles (200 each): histpp against both nulls; per-cell-type history against the distance-matched null.
#  6. Paired bootstrap (200): histpp, per-cell-type history.
#  7. Tables: notebooks/12_results_tables.md.
# Launch from the repo root:
#   nohup caffeinate -i scripts/run_session12.sh > data/session12.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }
H="TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist TRUNC_LEVELS=202304,202511,202602"
PP="TRUNC_INDEX=index_histpp.csv TRUNC_TAG=histpp"

stamp "1 missing vs attached synapses";   $PY scripts/missing_synapses.py
stamp "2 postsynaptic rollback: fetch";   $PY scripts/postsyn_history.py fetch
stamp "2 postsynaptic rollback: build";   $PY scripts/postsyn_history.py build
stamp "3 count-matched random loss sets"; $PY scripts/histsub_sets.py --reps 20

stamp "4 fits: histsub";  TRUNC_INDEX=index_histsub.csv TRUNC_TAG=histsub $R scripts/like_to_like_trunc.R fit poisson 4
stamp "4 fits: histpp";   env $PP $R scripts/like_to_like_trunc.R fit poisson 4
for CT in 4P 23P; do
  stamp "4 fits: history, $CT axons"; env $H TRUNC_PRE_CT=$CT $R scripts/like_to_like_trunc.R fit poisson 4
done

stamp "5 histpp: distance-matched shuffle"; env $PP TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 202304,202511,202602
stamp "5 histpp: cell-type shuffle";        env $PP $R scripts/like_to_like_trunc.R perm poisson 200 4 202304,202511,202602
for CT in 4P 23P; do
  stamp "5 history, $CT: distance-matched shuffle"
  env $H TRUNC_PRE_CT=$CT TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 202304,202511,202602
done

stamp "6 histpp: bootstrap"; env $PP $R scripts/like_to_like_trunc.R boot poisson 200 4
for CT in 4P 23P; do
  stamp "6 history, $CT: bootstrap"; env $H TRUNC_PRE_CT=$CT $R scripts/like_to_like_trunc.R boot poisson 200 4
done

stamp "7 tables"; $PY scripts/summarize_session12.py
stamp "done"
