#!/usr/bin/env bash
# Session 7: full truncation pipeline. Run from the repo root; every R step resumes if interrupted.
# Rough times on an 8 GB machine (from a 2-core test run; scale by your worker count):
#   truncation sets ~2 min; fit poisson ~8 min; fit tweedie ~15 min;
#   perm poisson (200 x 3 levels) ~1 h at 4 workers; boot poisson (200) ~3.5 h at 4 workers
#   (one bootstrap sample = 22 sets x 4 fits, ~4 min per worker);
#   perm/boot tweedie are ~8x slower per fit and use 2 workers (1.7 GB each): run them overnight.
set -euo pipefail
R=.r-env/bin/Rscript

.venv/bin/python scripts/truncate_like_to_like.py --reps 5 --block-um 25 --seed 7
$R scripts/like_to_like_trunc.R fit poisson 4
$R scripts/like_to_like_trunc.R fit tweedie 2
.venv/bin/python scripts/summarize_truncation.py          # preview: Wald CIs, no shuffle yet

$R scripts/like_to_like_trunc.R perm poisson 200 4 1,0.5,0.3
$R scripts/like_to_like_trunc.R boot poisson 200 4
.venv/bin/python scripts/summarize_truncation.py

# overnight (optional before the preprint figure is final):
# $R scripts/like_to_like_trunc.R perm tweedie 200 2 1,0.5,0.3
# $R scripts/like_to_like_trunc.R boot tweedie 100 2
