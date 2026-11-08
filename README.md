# wiring-worth

How does incomplete proofreading bias the wiring rules we infer from a connectome? This repository holds the
analysis behind the preprint *Incomplete proofreading biases wiring rules inferred from the MICrONS connectome
through spatial sampling* (W. Switzer, 2026; link to follow).

It uses the MICrONS dataset (mouse visual cortex: an EM connectome co-registered with two-photon imaging),
datastack `minnie65_public`, materialization **v1822** (pinned in `config.py`), accessed with `caveclient`.

**Main findings**
- On 193 co-registered excitatory neurons with now-complete axons, the like-to-like rule for receptive-field (RF)
  distance looked 10–12% steeper in those axons' partially proofread states (Apr 2023 – Feb 2026), read from the
  chunkedgraph edit history. The orientation rule looked weaker.
- The RF inflation is spatial sampling. A shuffled-wiring control that keeps each axon's target-distance
  distribution removes it: the preference beyond distance is the same in partial and complete states.
- Radial truncation of complete axons reproduces the RF inflation; random cable loss of the same amount does not.
- The orientation rule looked weaker because the 4P axons that carry it were the least complete: random loss of
  each axon's synapses down to its historical count reproduces the weakening exactly.
- Orientation like-to-like beyond distance is carried by 4P axons; 23P→23P shows none, matching Ding et al. (2025)
  for L2/3. This holds on the newest release (v1926), which nearly doubles the 23P→23P synapses.
- A second, structural rule is biased too: partial 4P axons lacked their ascending branches, which erased their
  preference for 23P over 4P targets. A distance-matched comparison does not correct this; the null has to match
  the geometry of the missing cable.

## Figures

| figure | script | inputs |
|---|---|---|
| Fig. 1 `figures/08_history.pdf` | `scripts/figure_history.py` | `data/history/`, `data/truncation/trunc_*_hist*.csv`, `trunc_perm_*` |
| Fig. 2 `figures/07_truncation_shift.pdf` | `scripts/figure_truncation_shift.py` | `data/truncation/trunc_{slopes,boot,perm}_poisson*.csv` |
| Fig. S1 `figures/14_celltype_history.pdf` | `scripts/figure_celltype_history.py` | `data/celltype_history/` |
| Baseline `figures/06_like_to_like.png` | `scripts/summarize_like_to_like.py` | `data/like_to_like/` |

Results notes for each session are in `notebooks/NN_results.md`; the preprint's numbers come from
`notebooks/06`–`15` (supplementary tables S1–S14 map to these files as listed in the preprint).

## Reproducing

**Environments.**
```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name microns
# R 4 + glmmTMB 1.1.14, pinned (conda/mamba):
mamba create -p .r-env --file r-environment.lock.txt   # osx-arm64; r-environment.yml for other platforms
```

**CAVE access** (one time): accept the MICrONS terms of service at
https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept, then run `.venv/bin/python scripts/save_token.py`.
Tokens are stored by `caveclient` outside the repository.

**Pipeline**, from the repo root. Downloaded tables go to `data/`, which is not committed. Every R step resumes
from its output CSV.

| step | command | what it does |
|---|---|---|
| 1 | `.venv/bin/python notebooks/01_setup_and_explore.py` | download the v1822 tables |
| 2 | `.venv/bin/python scripts/build_usable_set.py` | 694 co-registered, proofread excitatory neurons |
| 3 | `scripts/pull_synapses.py`, `fetch_skeletons.py`, `fetch_io_synapses.py`, `relabel_skeletons.py`, `compute_proximity.py` | synapses, skeletons, axon/dendrite split, Ld |
| 4 | `scripts/prepare_like_to_like.py`, then `.r-env/bin/Rscript scripts/like_to_like.R real` and `scripts/run_like_to_like_perms.sh` | baseline GLMMs + cell-type shuffle |
| 5 | `scripts/run_truncation_overnight.sh`, `run_truncation_boot.sh`, `run_session7_finish.sh`, `run_session7_sensitivity.sh` | simulated truncation, nulls, bootstrap |
| 6 | `scripts/edit_history.py fetch` / `build`, `run_edit_history.sh`, `run_history_ci.sh`, `run_distance_null.sh` | real proofreading history, distance-matched shuffle |
| 7 | `scripts/run_distance_sensitivity.sh` | bin width, within cell type |
| 8 | `scripts/run_session10.sh` | Feb 2026 state, shuffle top-ups, match-rule sensitivity, Ld models |
| 9 | `scripts/run_session11.sh` | column-cylinder and path-radius truncation; v1926 status of the 193 axons |
| 10 | `scripts/run_session12.sh` | missing vs. attached synapses, count-matched random loss, postsynaptic rollback |
| 11 | `scripts/run_session13.sh` | like-to-like by presynaptic × postsynaptic cell type |
| 12 | `.venv/bin/python scripts/celltype_connectivity_history.py` | target-type preference across the history (Table S13) |
| 13 | `scripts/run_v1926.sh` (sets `MICRONS_VERSION=1926` itself) | complete-axon results on v1926 (Table S14) |
| 14 | `scripts/figure_history.py`, `scripts/figure_truncation_shift.py`, `scripts/figure_celltype_history.py` | figures |

Runtimes on an 8 GB Mac are given at the top of each `run_*.sh`.

## Layout

- `notebooks/` — results notes and the setup script
- `scripts/` — analysis and run scripts
- `figures/` — output figures
- `data/` — downloaded and derived tables (not committed; derived tables are archived with each release)

## Citation and license

Code: MIT (see `LICENSE`). Cite via `CITATION.cff`. MICrONS data: cite the MICrONS Consortium (2025), *Nature*,
doi:10.1038/s41586-025-08790-w, and follow the dataset's terms of use.
