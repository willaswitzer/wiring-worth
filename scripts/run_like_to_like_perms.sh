#!/bin/sh
# Session 6 permutation schedule (8 GB machine). Resumable: rerunning continues each file.
# Poisson track (Nsyn + density) is cheap -> 6 workers; Tweedie track (Ld) needs ~1.7 GB/fit -> 2 workers.
set -e
cd "$(dirname "$0")/.."
R=.r-env/bin/Rscript
$R scripts/like_to_like.R perm main     poisson 1000 20 6
$R scripts/like_to_like.R perm main     tweedie 1000 20 2
$R scripts/like_to_like.R perm soma     poisson 1000  8 6
$R scripts/like_to_like.R perm soma     tweedie 1000  8 2
$R scripts/like_to_like.R perm ext_dend poisson 1000  8 6
$R scripts/like_to_like.R perm ext_dend tweedie 1000  8 2
