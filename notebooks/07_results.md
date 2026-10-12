# Session 7 results: like-to-like vs. simulated axon completeness

Final, at v1822 on the pinned environment (`.r-env`, glmmTMB 1.1.14):
- **Fits and shuffles:** `scripts/run_truncation_overnight.sh`, 2026-10-09 (00:35 → 05:21).
- **Bootstrap with rotated null reps:** `scripts/run_truncation_boot.sh`, 2026-10-09/10.
- **Calibrated truncation, within-cell-type and partial-vs-full analyses:** `scripts/run_session7_finish.sh`, 2026-10-10 (18:56 → 20:18).
- **Preprint figure:** `figures/07_truncation_shift.png` (`scripts/figure_truncation_shift.py`).
- All fits converged: 0 failures across fits, shuffles and bootstrap samples.
- The first overnight run (2026-10-08) is kept in `data/truncation/run1/`. It fitted V1→HVA (FF) at every level, which made those fits degenerate (see Design). The V1→V1 results of the two runs agree to within 0.01.

**Pipeline** (`scripts/run_truncation_overnight.sh`):
1. `truncate_like_to_like.py` → 78 truncation sets
2. `like_to_like_trunc.R fit` (Poisson, Tweedie)
3. `perm` (shuffled wiring, 200 per level at f = 1, 0.5, 0.3; Poisson and Tweedie)
4. `boot` (paired bootstrap over presynaptic neurons, 200; Poisson)
5. `summarize_truncation.py` → `07_truncation_table.md`, `figures/07_truncation.png`

## Design

- **Reference set:** the 193 fully extended axons, with the Session 6 pairs and main-variant models.
- **Level f:** the fraction of each axon's own cable kept, from 1.0 down to 0.3.
- **radial:** keep cable within the path length that keeps a fraction f; Ld and synapses follow the kept cable.
- **blocks (null):** the same cable loss, removed as random 25 µm cubes (5 reps).
- **synsub:** each axon's synapses are subsampled at random to radial's count, and Ld is left complete (5 reps).
- **Bootstrap:** sample b uses rep (b − 1) mod 5 of each null, so the null CIs include the nulls' own randomness. An earlier rep-0-only bootstrap is kept in `data/truncation/run2_boot_rep0/`; it overstated the synsub shift, because rep 0 is its most extreme rep.
- **Projection types:** only those that pass Ding's inclusion rule on the real data at *every* radial level are fitted.
  - RF: V1→V1 and HVA→HVA. Orientation: V1→V1.
  - V1→HVA (FF) is excluded: its synapses are distal, so 119 → 3 remain at f = 0.5, and 0 for orientation at f = 0.3.
  - A bootstrap sample that leaves a projection without synapses drops it from that fit.
- **Checks:** at f = 1, the recomputed Ld and Nsyn equal Session 6 exactly. With the reduced projection set, the V1→V1 slopes are within 0.001 of Session 6.

**Synapses and cable kept:**

| f | 0.9 | 0.7 | 0.5 | 0.3 |
|---|---|---|---|---|
| radial: synapses onto usable neurons kept | 0.94 | 0.82 | 0.65 | 0.42 |
| blocks: synapses kept | 0.90 | 0.70 | 0.49 | 0.30 |

Co-registered targets sit near the presynaptic somata. So at the same f, blocks lose *more* synapses than radial truncation, and the blocks null is conservative on sample size.

## Results (V1→V1; slope per 10° on the log scale; negative = like-to-like)

**Shift relative to full axons, slope(f) − slope(1).** Full-data estimate (nulls: mean over 5 reps) [95% paired bootstrap over presynaptic neurons, 200 samples]:

| model, similarity | radial f = 0.7 | radial 0.5 | radial 0.3 | blocks 0.3 | synsub 0.3 |
|---|---|---|---|---|---|
| Nsyn, RF | −0.08 [−0.13, −0.02] | −0.14 [−0.22, −0.06] | **−0.28 [−0.41, −0.14]** | +0.01 [−0.15, +0.16] | −0.09 [−0.23, +0.04] |
| Nsyn/Ld, RF | −0.04 [−0.11, +0.03] | −0.06 [−0.16, +0.04] | −0.14 [−0.30, +0.01] | −0.00 [−0.19, +0.18] | −0.11 [−0.28, +0.07] |
| Nsyn, orientation | −0.009 [−0.019, −0.001] | −0.014 [−0.031, −0.003] | −0.017 [−0.041, +0.001] | +0.001 [−0.033, +0.036] | +0.003 [−0.031, +0.034] |
| Nsyn/Ld, orientation | −0.008 [−0.017, 0.000] | −0.011 [−0.028, 0.000] | −0.013 [−0.037, +0.005] | +0.004 [−0.026, +0.035] | +0.002 [−0.032, +0.034] |

The bootstrap of the Ld (Tweedie) models hasn't run yet, so the Ld bands in the figure are Wald CIs.

**Paired contrasts** (full-data estimate [95% paired bootstrap]):

| contrast | f = 0.7 | f = 0.5 | f = 0.3 |
|---|---|---|---|
| radial − blocks, RF Nsyn | −0.08 [−0.20, +0.02] | −0.13 [−0.27, +0.01] | **−0.28 [−0.47, −0.11]** |
| radial − synsub, RF Nsyn | −0.05 [−0.13, +0.01] | −0.10 [−0.21, −0.00] | **−0.19 [−0.36, −0.02]** |
| radial − synsub, orientation Nsyn | −0.011 [−0.028, +0.004] | −0.017 [−0.041, +0.004] | −0.020 [−0.053, +0.015] |
| radial − synsub, orientation Nsyn/Ld | −0.010 [−0.027, +0.006] | −0.015 [−0.038, +0.006] | −0.015 [−0.048, +0.017] |

**Against the shuffled-wiring control** (radial; real slope vs. shuffled mean [95% interval], p; 200 shuffles each, so p = 0.005 is the floor):

| model, similarity | f = 1 | f = 0.5 | f = 0.3 |
|---|---|---|---|
| Ld, RF | −0.897 vs −0.166 [−0.26, −0.07], p = 0.005 | −1.019 vs −0.174, p = 0.005 | −1.104 vs −0.189 [−0.32, −0.04], p = 0.005 |
| Nsyn, RF | −0.856 vs −0.259 [−0.41, −0.12], p = 0.005 | −0.992 vs −0.295, p = 0.005 | −1.132 vs −0.329 [−0.52, −0.13], p = 0.005 |
| Nsyn/Ld, RF | −0.096 vs −0.092, p = 0.95 | −0.150 vs −0.124, p = 0.70 | −0.238 vs −0.136 [−0.29, +0.01], p = 0.20 |
| Ld, orientation | −0.018 vs −0.001, p = 0.005 | −0.019 vs 0.000, p = 0.010 | −0.021 vs 0.000, p = 0.010 |
| Nsyn, orientation | −0.040 vs +0.002, p = 0.005 | −0.055 vs +0.003, p = 0.005 | −0.057 vs +0.002, p = 0.005 |
| Nsyn/Ld, orientation | −0.021 vs +0.004, p = 0.055 | −0.032 vs +0.006, p = 0.020 | −0.034 vs +0.004, p = 0.050 |

HVA→HVA (RF): Ld goes −0.86 → −1.30 (p = 0.005 at every level), and Nsyn −0.57 → −0.90 (p = 0.080 → 0.040). Same direction as V1→V1, but only 89 → 38 synapses.

## Interpretation

1. **Radial truncation inflates the RF like-to-like rule; random cable loss doesn't.**
   - At 30% of cable, the RF synapse-count slope steepens by −0.28 [−0.41, −0.14], about 32%. Contact (Ld) follows: −0.90 → −1.10.
   - Random cable loss of the same amount doesn't shift it: +0.01 [−0.15, +0.16]. The radial − blocks contrast excludes 0: −0.28 [−0.47, −0.11].
   - The bias is already detectable at 70% of cable: −0.08 [−0.13, −0.02].
   - **But it comes from 4P axons, not 23P** (see the cell-type section). Real partial 23P axons sit at a cable ratio of ~0.75 (Session 5), yet 23P shows no detectable bias.
2. **Part of the inflation is shared with the shuffled null; most is not.**
   - The shuffled RF slope also steepens with truncation (Nsyn −0.26 → −0.33; Ld −0.17 → −0.19). Truncation strengthens the between-group RF–connectivity association that the shuffle keeps (Session 6, point 1).
   - The real-minus-shuffle gap still grows: Nsyn −0.60 → −0.80, Ld −0.73 → −0.92. That is about 75% (Nsyn) and 89% (Ld) of the raw shift.
   - The preprint should report both the raw shift and the change in the gap.
3. **Truncation doesn't create an RF rule in synapses per contact (Nsyn/Ld).**
   - The slope drifts more negative (−0.10 → −0.24), but stays inside the shuffle (p = 0.20 at f = 0.3), and its shift CI includes 0.
   - So the Session 6 conclusion holds under truncation: the RF rule is in contact, not in how contact is converted to synapses.
   - Caveat for partially proofread data: at 30% of cable the point estimate is 2.5× its full-axon value.
4. **Most of the radial effect is within-axon and spatial; little of it is re-weighting across axons.**
   - Synsub keeps each axon's synapse count as radial truncation does, but removes synapses at random. It shifts RF Nsyn by −0.09 [−0.23, +0.04], which can't be distinguished from 0.
   - The within-axon spatial part (radial − synsub) is −0.19 [−0.36, −0.02].
   - Cell-type composition doesn't explain the effect either. Radial truncation does keep more 23P/4P than 5P-IT synapses (48%, 41%, 27%), but re-weighting the full-axon per-cell-type slopes (23P −1.09, 4P −0.78, 5P-IT −0.63) by the synapses kept at f = 0.3 moves the pooled slope by only −0.02. The effect is *within* cell types (see the cell-type section).
   - For orientation, no component is separable from 0.
5. **Orientation:** the rule is beyond the shuffle at every level for Ld and Nsyn.
   - Truncation strengthens Nsyn modestly (−0.040 → −0.057). The shift CI excludes 0 at f = 0.7 and 0.5 and is borderline at 0.3; none of the contrasts with the nulls excludes 0.
   - Nsyn/Ld becomes nominally significant at f = 0.5 (p = 0.020), against p = 0.055 for full axons. A partially proofread dataset could therefore report an orientation rule at the synaptic scale that full axons only weakly support.
6. **Truncation erases interareal (feedforward) wiring first.** V1→HVA synapses are distal: 119 → 3 at f = 0.5. Rules for interareal projections can't be measured from radially truncated axons at all. HVA→HVA survives (89 → 38 synapses at f = 0.3) and shows the same RF inflation as V1→V1.

## Cell-type check: the bias is in 4P, not 23P

`scripts/celltype_truncation_check.R`: the V1→V1 RF Nsyn model is fitted separately to 23P and 4P presynaptic axons, with a paired bootstrap over that cell type's neurons (200 samples). Rerun on the pinned environment, which reproduced the earlier cloud-test numbers exactly.

| cell type | slope at f = 1 | shift at f = 0.7 | shift at f = 0.5 | shift at f = 0.3 |
|---|---|---|---|---|
| 23P (64 axons, 940 synapses) | −1.09 | −0.04 [−0.16, +0.03] | −0.05 [−0.22, +0.06] | −0.03 [−0.27, +0.14] |
| 4P (90 axons, 2,063 synapses) | −0.78 | −0.07 [−0.15, −0.01] | −0.12 [−0.22, −0.02] | **−0.33 [−0.50, −0.18]** |

5P-IT (21 axons) moves from −0.63 to −1.08 at f = 0.3, but its SE is 0.34, so it can't be interpreted.

- **The pooled radial effect comes from 4P (and possibly 5P-IT) axons.** 4P's RF preference is concentrated near the soma, and its distal synapses go to RF-dissimilar targets. 23P axons show the same RF preference near and far, so cutting them back barely changes their slope.
- **This matters for what the preprint can claim.** Session 5 found that real partially extended *23P* axons look radially truncated, while real partial *4P* axons look thinned beyond ~200 µm rather than cut off.
  - So the cell type whose real incompleteness matches the radial model (23P) shows no detectable bias. Its CI at f = 0.7 is [−0.16, +0.03], so a bias of up to about 15% can't be excluded.
  - The cell type that would be strongly biased by radial truncation (4P) is, in the real data, incomplete in a different way.
- **Revised framing:** "incompleteness biases the inferred rule" depends on how incompleteness interacts with each cell type's spatial wiring layout. Simulations need to match each cell type's real incompleteness pattern before predicting the bias in a real partially proofread dataset.
- This motivated the two analyses below: a **calibrated** truncation that matches real partial axons per cell type, and the plan's **third curve** (real partial vs. full axons).

## Calibrated truncation: realistic incompleteness barely biases the rule

`scripts/calibrated_truncation.py`: each fully extended 23P and 4P axon gets random **branch cuts**, and each cut removes everything beyond it.
- **Calibration:** the cut probability per path length is set so the expected kept cable at each distance equals the real partial/full cable ratio for that cell type (Session 5 profile).
- **Fit to target:** the achieved profile matches within 0.06 in every bin. Mean cable kept is 0.745 for 23P (real ratio ~0.75) and 0.78 for 4P (~0.80).
- **Other cell types:** no calibration exists, so their 26 axons are left complete.
- **Repeats and null:** 20 random repeats, plus a null (`calibblk`) that removes the same cable per axon as random 25 µm cubes.
- **Fitting and uncertainty:** `like_to_like_trunc.R` with `TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib`. The paired bootstrap rotates through the 20 repeats.
- **Environment:** pinned (`scripts/run_session7_finish.sh`). The numbers are identical to the earlier cloud test.

**Shift in the V1→V1 slope vs. full axons**, mean over repeats [95% paired bootstrap, 200 samples]:

| | full-axon slope | calibrated | cable-matched random null |
|---|---|---|---|
| RF, Nsyn (pooled) | −0.856 | **−0.012 [−0.090, +0.073]** | +0.011 [−0.072, +0.105] |
| RF, Nsyn/Ld | −0.096 | +0.009 [−0.062, +0.097] | +0.017 [−0.074, +0.131] |
| orientation, Nsyn | −0.040 | −0.001 [−0.016, +0.013] | −0.002 [−0.019, +0.014] |
| orientation, Nsyn/Ld | −0.021 | −0.001 [−0.016, +0.011] | −0.003 [−0.018, +0.012] |
| RF, Nsyn, 23P only | −1.091 | −0.015 (rep SD 0.04) | −0.013 |
| RF, Nsyn, 4P only | −0.783 | −0.021 (rep SD 0.05) | +0.006 |

- **The RF rule is about 1% steeper on average, and the CI rules out a bias beyond about ±10%.** For comparison, uniform radial truncation to 30% of cable shifts it by −0.27, about 32%.
- **Orientation:** the calibrated shift is near zero. Its CI (±0.015) is wide relative to the small orientation slope (−0.040).
- **Calibrated vs. null:** the difference is −0.025 [−0.109, +0.068] for RF Nsyn.
- **Shuffled-wiring control** (100 shuffles, reps rotated): under calibrated incompleteness the RF Nsyn slope is −0.868 vs. a shuffled mean of −0.266, p = 0.010 (the floor for 100 shuffles). For full axons it is −0.856 vs. −0.264. Orientation Nsyn: −0.041 vs. +0.003, p = 0.020.

## Real partial vs. full axons (the plan's third curve)

`scripts/prepare_partial.py` + `scripts/partial_vs_full.R`: the same models fitted to real partially extended axons (23P: 107; 4P: 167) and fully extended axons (23P: 64; 4P: 90), V1→V1.
- **Model:** `y ~ sim × axon_level + covariates + (1 | pre)`.
- **Uncertainty:** bootstrap over presynaptic neurons within cell type × axon level (200 samples).
- **Shuffled control:** 100 shuffles.

| | full (p vs shuffle) | partial (p vs shuffle) | partial − full [95% bootstrap] | p of difference vs shuffle | + soma distance |
|---|---|---|---|---|---|
| 23P, RF, Nsyn | −1.08 (0.010) | −1.14 (0.010) | −0.06 [−0.36, +0.26] | 0.91 | −0.05 |
| 4P, RF, Nsyn | −0.79 (0.010) | −0.93 (0.010) | −0.14 [−0.37, +0.15] | 0.52 | −0.04 |
| 23P, RF, Nsyn/Ld | −0.22 (1.00) | −0.26 (0.62) | −0.04 [−0.36, +0.30] | 0.76 | −0.09 |
| 4P, RF, Nsyn/Ld | −0.13 (0.49) | +0.05 (0.39) | +0.17 [−0.06, +0.42] | 0.14 | +0.09 |
| 23P, orientation, Nsyn | −0.00 (0.79) | −0.06 (0.020) | −0.054 [−0.119, +0.017] | 0.20 | −0.054 |
| 4P, orientation, Nsyn | −0.06 (0.010) | −0.07 (0.010) | −0.016 [−0.079, +0.042] | 0.86 | −0.020 |

p = 0.010 is the floor for 100 shuffles.

- **RF:** real partial and full axons show the same rule within error, consistent with the calibrated prediction (−0.015 for 23P, −0.021 for 4P).
  - **This is a weak test:** the CIs are about ±0.3, roughly ±30% of the rule. The real data can't rule out a bias of that size by itself; the calibrated simulation is what bounds it to about ±10%.
  - The 4P point estimate (−0.14) shrinks to −0.04 with soma distance. Partial and full 4P somata sit differently relative to the co-registered targets (Session 5): partial 4P axons have less than half the Ld per axon (9.5 vs 21.6 mm).
- **Orientation in 23P:** partial 23P axons show an orientation rule beyond the shuffle (p = 0.020), and full 23P axons don't (p = 0.79). But the difference is not significant (−0.054 [−0.119, +0.017], p = 0.20). Not a finding; at most a lead.

## Sensitivity (2026-10-11 overnight, `scripts/run_session7_sensitivity.sh`, pinned environment)

All fits converged.

**Calibration uncertainty.** `calibrated_truncation.py --profile-boot`: each of 20 reps draws its own target profile by resampling real partial and full 23P/4P axons (within cell type × axon level). Each draw is then used for that rep's branch cuts. The paired bootstrap rotates through the reps, so the CI includes uncertainty in the calibration, the random cuts and the sampling of presynaptic neurons.
- **Drawn profiles** (mean ± SD across draws, by path bin 0–200 … 800+ µm):
  - 23P: 0.92 ± 0.03, 0.82 ± 0.04, 0.57 ± 0.04, 0.54 ± 0.05, 0.42 ± 0.09.
  - 4P: 0.93 ± 0.03, 0.75 ± 0.04, 0.65 ± 0.05, 0.74 ± 0.06, 0.73 ± 0.23.
- **RF Nsyn shift:** −0.016 [−0.100, +0.046]. The fixed profile gave −0.012 [−0.090, +0.073], and the SD across reps barely changes (0.029 vs 0.030), so the calibration's own uncertainty adds almost nothing. The bound on the realistic bias stays at about ±12% of the rule.
- **Other models:** RF Nsyn/Ld +0.004 [−0.079, +0.082]; orientation Nsyn +0.001 [−0.015, +0.016].
- The two test draws earlier (shift ≈ −0.055) were at the tail of this distribution.

**Block size of the random-cable-loss null.** 10 µm and 50 µm cubes, against 25 µm in the main run:

| null cube size | RF Nsyn shift at f = 0.3 | radial − blocks at f = 0.3 (paired) | radial − blocks at f = 0.5 |
|---|---|---|---|
| 10 µm | −0.003 [−0.24, +0.19] | −0.27 [−0.52, −0.02] | −0.13 [−0.27, +0.04] |
| 25 µm (main) | +0.01 [−0.15, +0.16] | −0.28 [−0.47, −0.11] | −0.13 [−0.27, +0.01] |
| 50 µm | −0.023 [−0.21, +0.18] | −0.25 [−0.47, −0.01] | −0.11 [−0.26, +0.04] |

- **The null stays flat at every cube size and level** (|shift| ≤ 0.03 at every f).
- **Radial truncation differs from it at f = 0.3 at every cube size.** The CIs are wider for 10 and 50 µm, because those nulls vary more between reps.
- **Pairing check:** the bootstraps share presynaptic draws across runs (same seeds), and the f = 1 bootstrap slopes are identical across files, which is what makes the radial − blocks contrast across runs valid.
- **Orientation:** no contrast excludes 0.

## Revised summary for the preprint

1. Incompleteness *can* bias an inferred wiring rule, but how much depends on the pattern of loss and on each cell type's spatial wiring layout.
   - Cutting axons back uniformly from the far end inflates the RF like-to-like rule by up to 32%. The effect comes from 4P axons, whose RF-similar partners lie near the soma.
   - Random loss of the same amount of cable has no effect.
2. Simulated with MICrONS's actual partial-proofreading profiles, the bias is ~1-2%: −0.016 [−0.100, +0.046], including uncertainty in the calibration itself. That is within about ±12% of the rule.
3. Real partial axons agree with full ones, as the calibrated simulation predicts, but within wide error (about ±30%). The calibrated simulation is what bounds the bias to about ±10%.
4. **So the published like-to-like rule is robust to the reconstruction incompleteness present in MICrONS.**
   - The method (calibrate a cut model on partial vs. full axons, then simulate) gives an error budget for any wiring rule, which mirrors the Fly project.
   - Interareal (V1→HVA) projections are the exception: they are lost first under cutting back from the far end.

## Next

- **Preprint figure:** `figures/07_truncation_shift.png`.
  - Columns 1–2: shift vs. completeness for radial truncation and both nulls, with the shift beyond the shuffle.
  - Column 3: calibrated prediction vs. real partial − full axons.
- **Optional:** Ld (Tweedie) bootstrap. Calibration uncertainty and block size are done (see Sensitivity).
- **Write the preprint around the revised summary above.**
