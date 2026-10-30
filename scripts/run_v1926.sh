#!/bin/sh
# v1926 add-on: complete-axon like-to-like on the newest public release, next to the pinned v1822 results.
# v1926 adds up to 49 co-registered 23P neurons with fully extended axons (before the quality filter), roughly
# doubling the 23P axons; the question is whether 23P->23P orientation stays flat with more data.
# Rebuilds the usable set at v1926 and reruns only the complete-axon analyses: baseline + cell-type shuffle, the
# distance-matched shuffle, per presynaptic cell type, and per pre x post cell type. No truncation or edit history.
# All outputs carry a _v1926 suffix (or live in data/*_v1926/); the pinned v1822 outputs are not touched.
# Downloads tables, synapses and skeletons for new roots (skeleton cache is shared, keyed by root id).
# Every step resumes if rerun. Launch from the repo root:
#   nohup caffeinate -i scripts/run_v1926.sh > data/v1926.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }
export MICRONS_VERSION=1926

stamp "tables";              $PY notebooks/01_setup_and_explore.py
stamp "usable set";          $PY scripts/build_usable_set.py
stamp "synapses";            $PY scripts/pull_synapses.py
stamp "skeletons";           $PY scripts/fetch_skeletons.py
stamp "skeletons (retry)";   $PY scripts/fetch_skeletons.py
stamp "axon/dendrite synapses"; $PY scripts/fetch_io_synapses.py || $PY scripts/fetch_io_synapses.py
stamp "relabel new roots";   $PY scripts/relabel_skeletons.py
stamp "proximity";           $PY scripts/compute_proximity.py
stamp "pairs";               $PY scripts/prepare_like_to_like.py
stamp "complete-axon set";   $PY scripts/make_full_set.py

stamp "baseline fits";       $R scripts/like_to_like.R real
stamp "cell-type shuffle";   $R scripts/like_to_like.R perm main poisson 200 600 4

export TRUNC_LEVELS=1
stamp "pooled: distance-matched"; TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 1
for PRE in 4P 23P; do
  export TRUNC_PRE_CT=$PRE
  stamp "$PRE axons: fit";              $R scripts/like_to_like_trunc.R fit poisson 4
  stamp "$PRE axons: distance-matched"; TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 1
done
for C in 23P:23P 4P:4P 4P:23P 23P:4P; do
  PRE=${C%%:*}; POST=${C##*:}
  export TRUNC_PRE_CT=$PRE TRUNC_POST_CT=$POST
  stamp "$PRE -> $POST: fit";              $R scripts/like_to_like_trunc.R fit poisson 4
  stamp "$PRE -> $POST: distance-matched"; TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 1
  stamp "$PRE -> $POST: cell-type shuffle"; $R scripts/like_to_like_trunc.R perm poisson 200 4 1
done
unset TRUNC_PRE_CT TRUNC_POST_CT TRUNC_LEVELS MICRONS_VERSION
stamp "table";               $PY scripts/summarize_v1926.py
stamp "done"
