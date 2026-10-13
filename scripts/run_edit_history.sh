#!/bin/sh
# Session 8: edit-history natural experiment (needs CAVE access; your Mac).
#   1. fetch: supervoxels of the 193 axons' output synapses + roots at monthly time points 2021-01 .. v1822
#      + axon labels at each queryable version (resumable; roughly 15-60 min, depending on the CAVE server)
#   2. build: each axon's state at each version -> truncation-format sets (~2 min)
#   3. fit the rule on each time point's state (~10 min). The bootstrap is run afterwards on a chosen
#      subset of time points (a bootstrap over all ~67 monthly points would take ~9 h).
# Quick test first (2 axons, 2 versions, writes data/history_test/ only):
#   .venv/bin/python scripts/edit_history.py fetch --limit 2 --versions 1078,1822
# Launch:  nohup caffeinate -i scripts/run_edit_history.sh > data/edit_history.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
PY=.venv/bin/python
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }
stamp "fetch";  $PY scripts/edit_history.py fetch
stamp "build";  $PY scripts/edit_history.py build
stamp "fit";    TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist $R scripts/like_to_like_trunc.R fit poisson 4
stamp "done"
