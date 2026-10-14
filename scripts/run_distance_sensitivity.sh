#!/bin/sh
# Session 9: sensitivity of the distance-matched shuffle (open items in the preprint).
#  1. Bin width: 25 and 100 um soma-distance bins (50 um is the main analysis), on complete axons,
#     radial f = 0.3 and the Nov 2025 history state. Outputs get a _bin25 / _bin100 suffix.
#  2. Within cell type (50 um bins): 4P and 23P presynaptic axons separately, complete and f = 0.3,
#     to test whether the 4P radial inflation is absorbed by distance within 4P itself.
# ~2.5 h at 4 workers. Launch:
#   nohup caffeinate -i scripts/run_distance_sensitivity.sh > data/distance_sensitivity.log 2>&1 &
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
stamp() { echo "=== $(date '+%H:%M:%S') $*"; }

for B in 25 100; do
  stamp "bin $B um: radial";  TRUNC_NULL=distance TRUNC_DIST_BIN_UM=$B $R scripts/like_to_like_trunc.R perm poisson 200 4 1,0.3
  stamp "bin $B um: history"; TRUNC_NULL=distance TRUNC_DIST_BIN_UM=$B TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist \
                                TRUNC_LEVELS=202511,1822 $R scripts/like_to_like_trunc.R perm poisson 200 4 202511
done

for CT in 4P 23P; do
  stamp "$CT: real fits";        TRUNC_PRE_CT=$CT TRUNC_LEVELS=1,0.3 $R scripts/like_to_like_trunc.R fit poisson 4
  stamp "$CT: distance null";    TRUNC_PRE_CT=$CT TRUNC_NULL=distance $R scripts/like_to_like_trunc.R perm poisson 200 4 1,0.3
done
stamp "done"
