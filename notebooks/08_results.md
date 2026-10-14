# Session 8 results: real proofreading history and a distance-matched null

Run 2026-10-13 at v1822 on the pinned environment (`.r-env`, glmmTMB 1.1.14). Both changes come from C. Schneider-Mizell's reply (Allen Institute, 2026-10-13):
- **Edit history:** partial axon extension in MICrONS followed axons to a fixed radial distance around the proofread column, so the bias of real incompleteness should be measured on the same cells' own past states, read from the chunkedgraph history.
- **Distance-matched null:** the shuffled-wiring control should keep the targets' spatial distribution, because nearby cells share thalamic input and receptive fields.

All fits converged: 68 history states, 200 bootstrap samples, 200 shuffles per condition and level.

**Scripts:**
- `scripts/edit_history.py fetch` / `build`
- `scripts/run_edit_history.sh`
- `scripts/run_distance_null.sh`
- `scripts/run_history_ci.sh`
- `like_to_like_trunc.R` options: `TRUNC_NULL=distance`, `TRUNC_LEVELS`

## 1. Real proofreading history of the 193 fully extended axons

For the first of every month from 2021-01 to v1822, each output synapse of each axon was counted as "attached" if its presynaptic supervoxel had the same root as the neuron's nucleus supervoxel at that time (`chunkedgraph.get_roots`). Each axon was then cut down to its attached cable (a 1 µm piece is kept if its nearest output synapse is attached). Ld and Nsyn were recomputed, and the Session 6 main-variant models were refitted.
- **Check:** the v1822 state reproduces Session 6 exactly.
- **Full table:** `notebooks/08_history_summary.md`.

| period | unproofread (< 10% of cable) | partial | complete (> 95%) | median cable attached | V1→V1 synapses attached (of 3,586) |
|---|---|---|---|---|---|
| Jan 2021 | 134 | 56 | 3 | 4% | 277 |
| Mar 2023 (proofreading push) | 25 | 147 | 21 | 33% | 1,758 |
| Apr 2023 – Nov 2025 (plateau) | 11–16 | 136–157 | 20–46 | 36–50% | 1,934–2,331 |
| Feb 2026 (completion push) | 4 | 65 | 124 | 100% | 2,744 |
| Jun 2026 (v1822) | 0 | 0 | 193 | 100% | 3,586 |

- **The plateau is the partial state Casey described.** About 150 axons sat at roughly 40% of their final cable for more than two years. They were completed between December 2025 and June 2026.
- **Labels confirm it.** Axon labels could still be queried at v943, v1300, v1507, v1621, v1718 and v1822. The number labelled `axon_fully_extended` rose from 40 at v1300 to 77 at v1621, 146 at v1718 and 193 at v1822.

**Shift of the V1→V1 slope vs. the complete state.** Each value is full-data minus complete-state, with a 95% paired bootstrap over presynaptic neurons (200 samples):

| state | RF, Nsyn | RF, Nsyn/Ld | orientation, Nsyn | orientation, Nsyn/Ld |
|---|---|---|---|---|
| Apr 2023 (plateau start) | −0.093 [−0.197, +0.002] | −0.022 [−0.149, +0.100] | +0.004 [−0.021, +0.023] | +0.007 [−0.016, +0.026] |
| Nov 2025 (plateau end) | **−0.106 [−0.197, −0.019]** | −0.054 [−0.150, +0.042] | +0.012 [−0.004, +0.026] | +0.012 [−0.003, +0.025] |
| Feb 2026 (mid-completion) | **−0.087 [−0.155, −0.021]** | −0.037 [−0.103, +0.042] | **+0.017 [+0.003, +0.030]** | **+0.016 [+0.003, +0.028]** |

Complete state (v1822): RF Nsyn −0.856, RF Nsyn/Ld −0.096, orientation Nsyn −0.040, orientation Nsyn/Ld −0.021.

- **RF:** while the axons were partially proofread, the RF synapse rule looked 10–12% steeper than on the finished axons. The CI excludes 0 at the plateau's end and during completion, and is borderline at its start.
- **Orientation:** in the partial-era data the orientation rule looked weaker, by up to 40% during completion (−0.023 vs −0.040).

## 2. Distance-matched shuffled-wiring control

The old null permuted postsynaptic functional properties among neurons of the same cell type × area. The new one permutes them, separately for each presynaptic neuron, only among candidate targets of the same cell type × area **and the same 50 µm soma-distance bin**. The shuffled targets therefore keep the real targets' spatial distribution.

**V1→V1, real slope vs. shuffled mean [95%], p (200 shuffles; floor 0.005):**

| | distance-matched null | cell-type null (Session 6/7) |
|---|---|---|
| RF Nsyn, complete axons | −0.856 vs −0.716 [−0.806, −0.621], p = 0.005 | −0.856 vs −0.259, p = 0.005 |
| RF Nsyn, radial f = 0.5 | −0.992 vs −0.906, p = 0.12 | −0.992 vs −0.295, p = 0.005 |
| RF Nsyn, radial f = 0.3 | −1.132 vs −1.010, p = 0.085 | −1.132 vs −0.329, p = 0.005 |
| RF Nsyn, history Apr 2023 | −0.949 vs −0.838, p = 0.080 | −0.949 vs −0.301, p = 0.005 |
| RF Nsyn, history Nov 2025 | −0.962 vs −0.814, p = 0.025 | −0.962 vs −0.308, p = 0.005 |
| RF Nsyn, calibrated cuts | −0.868 vs −0.738, p = 0.025 | (Session 7) |
| RF Nsyn/Ld, complete axons | −0.096 vs −0.007, p = 0.12 | −0.096 vs −0.092, p = 0.95 |
| Orientation Nsyn, complete axons | −0.040 vs −0.004, p = 0.005 | −0.040 vs +0.002, p = 0.005 |
| Orientation Nsyn, radial f = 0.3 | −0.057 vs −0.020, p = 0.035 | −0.057 vs +0.002, p = 0.005 |
| Orientation Nsyn, history Apr 2023 | −0.037 vs −0.013, p = 0.11 | −0.037 vs +0.001, p = 0.055 |
| Orientation Nsyn, history Nov 2025 | −0.028 vs −0.005, p = 0.10 | −0.028 vs +0.005, p = 0.035 |
| Orientation Nsyn/Ld, complete axons | −0.021 vs +0.008, p = 0.020 | −0.021 vs +0.004, p = 0.055 |

## Interpretation

1. **Most of the RF like-to-like rule is distance.** Nearby neurons have nearby receptive fields (retinotopy). Matching shuffled targets on soma distance moves the null from −0.26 to −0.72, so about 80% of the RF synapse slope is explained by where targets sit. The remaining like-to-like (−0.14) is still beyond every shuffle (p = 0.005).
2. **Incompleteness inflates the RF rule only through that distance effect.** The real-minus-distance-null gap is the same in every state:

    | state | real − distance-matched null |
    |---|---|
    | complete axons | −0.140 |
    | history, Apr 2023 | −0.112 |
    | history, Nov 2025 | −0.148 |
    | calibrated cuts | −0.130 |
    | radial, f = 0.5 | −0.086 |
    | radial, f = 0.3 | −0.122 |

    Partial axons over-sample nearby targets, so the raw slope steepens (by about 10% in the real history, up to 32% under strong radial truncation). The RF wiring preference beyond distance does not change.
3. **Incompleteness weakens the orientation rule and costs power.** Orientation like-to-like is beyond the distance-matched null on complete axons (p = 0.005). On the plateau-era axons it is weaker (−0.028 to −0.037) and not significant against either null (p = 0.10–0.11 distance-matched; 0.035–0.055 cell-type). A study run on the 2023–2025 data would have reported at most a borderline orientation rule.
4. **What this replaces in Session 7:**
    - The ~2% "realistic bias" from random branch cuts was a model of the wrong mechanism. The real history gives about 10% for RF.
    - Pooled radial truncation's inflation is absorbed by the distance-matched null, consistent with the 4P/thalamic-input explanation. A distance-matched null within 4P and 23P separately has not been run.
5. **For the preprint:**
    - The headline becomes "proofreading incompleteness inflated the RF rule by ~10% through spatial sampling and weakened the orientation rule; a distance-matched null removes the first and is needed to see the second."
    - The baseline must be reported against the distance-matched null as well as the cell-type null.

## Figure

`figures/08_history.png` (`scripts/figure_history.py`; values in `notebooks/08_figure_table.md`). It is Figure 1 of the preprint; the Session 7 truncation figure becomes Figure 2.

- Radial truncation strengthens the orientation rule (gap −0.048 at f = 0.5), while the real history weakened it, so a radial cutoff does not explain the orientation effect.

## Caveats

- **Only the presynaptic side is rolled back.** Synapse attribution on the postsynaptic side is today's, and each axon piece takes the state of its nearest output synapse.
- **Bin width:** the distance-matched null uses 50 µm bins; sensitivity to the bin width (25 and 100 µm) has not been run.
- **Bootstrap scope:** the history bootstrap covers four time points only (Apr 2023, Nov 2025, Feb 2026, v1822).
