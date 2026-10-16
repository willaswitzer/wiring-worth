# Sessions 10–11 results: remaining preprint analyses and C. Schneider-Mizell's follow-ups

Run 2026-10-15/16 at v1822 on the pinned environment. Scripts: `scripts/run_session10.sh` (12:54 → 21:49) and
`scripts/run_session11.sh` (13:09 → 00:05). All fits, shuffles and bootstrap samples converged. V1→V1, slope per
10°; p two-sided, centred on the shuffled mean.

## 1. Baseline with 500 shuffles (Table S1)

Poisson models topped up to 500 shuffles in every variant (Ld already had 500), so every p floor is now 0.002.

| variant | RF Nsyn (p) | RF Nsyn/Ld (p) | orientation Nsyn (p) | orientation Nsyn/Ld (p) |
|---|---|---|---|---|
| main | −0.855 vs −0.255 (0.002) | −0.096 vs −0.097 (1.00) | −0.040 vs +0.002 (0.002) | −0.021 vs +0.003 (0.084) |
| + soma distance | −0.642 vs −0.244 (0.002) | −0.167 vs −0.088 (0.21) | −0.036 vs +0.001 (0.014) | −0.026 vs +0.003 (0.042) |
| extended dendrites only | −0.970 vs −0.296 (0.002) | −0.102 vs −0.105 (0.96) | −0.047 vs +0.002 (0.004) | −0.022 vs +0.004 (0.080) |

The main-variant orientation Nsyn/Ld p moved from 0.061 (228 shuffles) to 0.084.

## 2. Co-registration match rules (Table S3)

`scripts/match_sensitivity.py` re-chooses each neuron's digital-twin unit within the main usable set; the main rule
reproduces the usable set exactly. Counts in `notebooks/10_match_rules.md`. Main variant, 200 shuffles:

| rule | neurons (fully extended pre) | RF Nsyn (p) | RF Nsyn/Ld (p) | orientation Nsyn (p) | orientation Nsyn/Ld (p) |
|---|---|---|---|---|---|
| main (score ≥ 0, max cc_max) | 694 (193) | −0.855 (0.002) | −0.096 (1.00) | −0.040 (0.002) | −0.021 (0.084) |
| lowest residual | 660 (186) | −0.862 (0.005) | −0.100 (0.995) | −0.051 (0.005) | −0.029 (0.015) |
| Ding match rule | 693 (193) | −0.851 (0.005) | −0.077 (0.76) | −0.037 (0.015) | −0.017 (0.10) |

The baseline does not depend on the match rule. 87 neurons would enter only under Ding's rule; they have no
proximity data and were not added.

## 3. February 2026 history state against both shuffles

| | real | distance-matched: gap (p) | cell-type: gap (p) |
|---|---|---|---|
| RF Nsyn | −0.943 | −0.156 (0.010) | −0.630 (0.005) |
| orientation Nsyn | −0.023 | −0.022 (0.10) | −0.028 (0.055) |

Same pattern as the plateau: RF beyond distance equals complete axons (−0.140); orientation weaker and not
distinguishable from the distance-matched shuffle. Figure 1 D/E regenerated with this state.

## 4. Ld (Tweedie) models (Table S6)

**Against the distance-matched shuffle** (200 shuffles):

| state | RF real | RF gap (p) | orientation real | orientation gap (p) |
|---|---|---|---|---|
| complete | −0.897 | −0.096 (0.005) | −0.018 | −0.007 (0.005) |
| radial f = 0.3 | −1.104 | −0.104 (0.005) | −0.021 | −0.005 (0.010) |
| Nov 2025 | −0.972 | −0.114 (0.005) | −0.019 | −0.006 (0.005) |

- About 90% of the RF contact rule and about 60% of the orientation contact rule are target distance.
- What remains beyond distance does not change with truncation or real incompleteness.
- So the orientation preference beyond distance is mainly in synapse count, not contact (Nsyn: −0.040 real,
  −0.037 beyond distance).

**Truncation bootstrap** (200 samples, paired over presynaptic neurons; replaces the Wald bands in
`figures/07_truncation.png`):

| | radial f = 0.5 | radial f = 0.3 | random cable loss f = 0.3 |
|---|---|---|---|
| RF, Ld shift | −0.122 [−0.151, −0.098] | −0.207 [−0.251, −0.175] | +0.044 [−0.002, +0.056] |
| orientation, Ld shift | −0.001 [−0.003, +0.001] | −0.003 [−0.006, +0.001] | −0.002 [−0.007, +0.002] |

Synapse subsampling leaves Ld complete by construction, so its Ld shift is 0.

## 5. Truncation at one radius for every axon (Session 11; Table S8)

Partial extension followed axons to a fixed radius around the proofread column (C. Schneider-Mizell). Full tables:
`notebooks/11_absolute_truncation.md`; set summary: `data/truncation/abs_summary_v1822.md`.

- **Column cylinder** (axis = median soma x, z of the 193 axons; 28 axons with somata > 100 µm away left complete,
  14 in V1). At Rc = 100 µm it keeps a median 38% of V1 cable and 68% of synapses, close to the real plateau
  (36–50% of cable, 54–65% of synapses).
  - RF Nsyn shift −0.174 [−0.272, −0.084]; beyond the distance-matched shuffle −0.149 (p = 0.010), against −0.140 for
    complete axons.
  - Orientation Nsyn shift −0.007 [−0.019, +0.008]; beyond distance −0.036 (p = 0.020).
- **Common path radius from the soma.** At 400 µm: RF shift −0.098 [−0.168, −0.038], beyond distance −0.111
  (p = 0.035). At 250 µm (41% of synapses) the RF shift is −0.251 and the gap beyond distance falls to −0.049
  (p = 0.50): power is lost. Orientation strengthens (−0.040 [−0.068, −0.018] at 250 µm).

**Conclusion.** Cutting at the column boundary, as proofreading did, inflates the raw RF rule (here by 20%; the
real history by 10–12%) and leaves the preference beyond distance unchanged. No geometric truncation weakens the
orientation rule, so the weaker orientation rule in the real partial states is not a geometric effect.

## 6. Newest release (Session 11; Table S9)

`scripts/check_new_release.py` → `notebooks/11_new_release.md`. The newest version is v1926 (materialized
2026-09-28). All 193 axons stay fully extended. It adds 2 co-registered 4P neurons with fully extended axons, and up to
49 23P and 21 5P-ET (before quality filtering and the dendrite requirement). The 4P results need no rerun; a v2 could
add the 23P axons.

## Figures

- `figures/08_history.pdf` (Figure 1): D/E now show complete, Apr 2023, Nov 2025, Feb 2026, the Rc = 100 µm column
  cylinder and radial f = 0.3. The calibrated branch cuts moved to the supplement.
- `figures/07_truncation_shift.pdf` (Figure 2): unchanged from 2026-10-15.
