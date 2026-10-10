# Session 7 results: like-to-like vs. simulated axon completeness

Final run 2026-10-09 at v1822 on the pinned environment (`.r-env`, glmmTMB 1.1.14), via `scripts/run_truncation_overnight.sh` (00:35 → 05:21).
- All fits converged: 0 failures across fits, 1,200 shuffles per track and 200 bootstrap samples.
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
- **Bootstrap:** uses rep 0 of each null.
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

**Shift relative to full axons, slope(f) − slope(1).** Paired bootstrap over presynaptic neurons, 200 samples, median [95% CI]:

| model, similarity | radial f = 0.7 | radial 0.5 | radial 0.3 | blocks 0.3 | synsub 0.3 |
|---|---|---|---|---|---|
| Nsyn, RF | −0.08 [−0.13, −0.02] | −0.13 [−0.22, −0.06] | **−0.27 [−0.41, −0.14]** | −0.05 [−0.19, +0.09] | −0.15 [−0.25, −0.02] |
| Nsyn/Ld, RF | −0.04 [−0.11, +0.03] | −0.06 [−0.16, +0.04] | −0.13 [−0.30, +0.01] | −0.10 [−0.29, +0.07] | −0.18 [−0.30, −0.01] |
| Nsyn, orientation | −0.010 [−0.019, −0.001] | −0.017 [−0.031, −0.003] | −0.019 [−0.041, +0.001] | −0.012 [−0.034, +0.018] | +0.018 [−0.002, +0.041] |
| Nsyn/Ld, orientation | −0.009 [−0.017, 0.000] | −0.014 [−0.028, 0.000] | −0.014 [−0.037, +0.005] | −0.008 [−0.029, +0.020] | +0.018 [−0.001, +0.041] |

The bootstrap of the Ld (Tweedie) models hasn't run yet, so the Ld bands in the figure are Wald CIs.

**Paired contrasts** (bootstrap median [95% CI]):

| contrast | f = 0.7 | f = 0.5 | f = 0.3 |
|---|---|---|---|
| radial − blocks, RF Nsyn | −0.08 [−0.16, 0.00] | −0.09 [−0.21, +0.03] | **−0.22 [−0.39, −0.05]** |
| radial − synsub, RF Nsyn | −0.03 [−0.10, +0.03] | −0.08 [−0.17, +0.01] | −0.13 [−0.27, +0.02] |
| radial − synsub, orientation Nsyn | **−0.019 [−0.029, −0.010]** | **−0.031 [−0.045, −0.015]** | **−0.037 [−0.061, −0.014]** |
| radial − synsub, orientation Nsyn/Ld | **−0.019 [−0.028, −0.008]** | **−0.028 [−0.041, −0.011]** | **−0.033 [−0.056, −0.009]** |

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
   - At 30% of cable, the RF synapse-count slope steepens by −0.27 [−0.41, −0.14], about 32%. Contact (Ld) follows: −0.90 → −1.10.
   - Random cable loss of the same amount doesn't shift it: −0.05 [−0.19, +0.09]. The radial − blocks contrast excludes 0: −0.22 [−0.39, −0.05].
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
4. **The synsub null drifts somewhat, but cell-type composition does not explain the radial effect.**
   - Subsampling each axon's synapses to radial's count steepens RF Nsyn by −0.09 at f = 0.3 (mean of 5 reps; reps range from −0.06 to −0.15).
   - The bootstrap CIs for both nulls above used rep 0 only, and for synsub rep 0 is the most extreme rep. They are being recomputed with the reps rotated across bootstrap samples (`scripts/run_truncation_boot.sh`). Radial CIs are unaffected.
   - **Composition is not the explanation.** Radial truncation does keep more 23P/4P than 5P-IT synapses (48%, 41%, 27%). But re-weighting the full-axon per-cell-type slopes (23P −1.09, 4P −0.78, 5P-IT −0.63) by the synapses kept at f = 0.3 moves the pooled slope by only −0.02. The radial effect is *within* cell types (see the cell-type section below).
5. **Orientation:** the rule is beyond the shuffle at every level for Ld and Nsyn.
   - Truncation strengthens Nsyn modestly (−0.040 → −0.057). The shift CI excludes 0 at f = 0.7 and 0.5 and is borderline at 0.3.
   - Nsyn/Ld becomes nominally significant at f = 0.5 (p = 0.020), against p = 0.055 for full axons. A partially proofread dataset could therefore report an orientation rule at the synaptic scale that full axons only weakly support.
6. **Truncation erases interareal (feedforward) wiring first.** V1→HVA synapses are distal: 119 → 3 at f = 0.5. Rules for interareal projections can't be measured from radially truncated axons at all. HVA→HVA survives (89 → 38 synapses at f = 0.3) and shows the same RF inflation as V1→V1.

## Cell-type check: the bias is in 4P, not 23P

`scripts/celltype_truncation_check.R`: the V1→V1 RF Nsyn model is fitted separately to 23P and 4P presynaptic axons, with a paired bootstrap over that cell type's neurons (200 samples). This was run in the cloud test environment (glmmTMB 1.1.8, which reproduced the pinned fits exactly at f = 1).

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
- **Next analyses:**
  - A **calibrated** condition: remove cable per cell type with the probability-vs-path-distance profile of real partial axons (Session 5 cable-deficit table, 23P and 4P). This predicts the bias in a dataset like the real one.
  - The plan's **third curve**: real 23P (and 4P) partial vs. full axons, as a direct test of that prediction.

## Calibrated truncation: realistic incompleteness barely biases the rule

`scripts/calibrated_truncation.py`: each fully extended 23P and 4P axon gets random **branch cuts**, and each cut removes everything beyond it.
- **Calibration:** the cut probability per path length is set so the expected kept cable at each distance equals the real partial/full cable ratio for that cell type (Session 5 profile).
- **Fit to target:** the achieved profile matches within 0.06 in every bin. Mean cable kept is 0.745 for 23P (real ratio ~0.75) and 0.78 for 4P (~0.80).
- **Other cell types:** no calibration exists, so their 26 axons are left complete.
- **Repeats and null:** 20 random repeats, plus a null (`calibblk`) that removes the same cable per axon as random 25 µm cubes.
- **Fitting and uncertainty:** `like_to_like_trunc.R` with `TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib`. The paired bootstrap rotates through the 20 repeats.
- **Environment:** run here in the cloud test environment (glmmTMB 1.1.8), to be repeated on the pinned environment by `scripts/run_session7_finish.sh`.

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

## Real partial vs. full axons (the plan's third curve)

`scripts/prepare_partial.py` + `scripts/partial_vs_full.R`: the same models fitted to real partially extended axons (23P: 107; 4P: 167) and fully extended axons (23P: 64; 4P: 90), V1→V1. The model is `y ~ sim × axon_level + covariates + (1 | pre)`. Wald estimates are below; the bootstrap and shuffled control are in the final run.

| | full | partial | partial − full (SE) | + soma distance |
|---|---|---|---|---|
| 23P, RF, Nsyn | −1.08 | −1.14 | −0.06 (0.14) | −0.05 (0.14) |
| 4P, RF, Nsyn | −0.79 | −0.93 | −0.14 (0.10) | −0.04 (0.10) |
| 23P, RF, Nsyn/Ld | −0.22 | −0.26 | −0.04 (0.16) | −0.09 (0.15) |
| 4P, RF, Nsyn/Ld | −0.13 | +0.05 | +0.17 (0.11) | +0.09 (0.11) |
| 23P, orientation, Nsyn | −0.00 | −0.06 | −0.054 (0.028) | −0.054 (0.028) |
| 4P, orientation, Nsyn | −0.06 | −0.07 | −0.016 (0.021) | −0.020 (0.021) |

- **RF:** real partial axons show the same rule as full axons within error. The calibrated simulation predicted shifts of −0.015 (23P) and −0.021 (4P), and the observed differences are consistent with both that prediction and with 0.
  - The 4P difference (−0.14) shrinks to −0.04 once soma distance is included. Partial and full 4P somata sit differently relative to the co-registered targets (Session 5): partial 4P axons have less than half the Ld per axon (9.5 vs 21.6 mm).
- **Orientation in 23P is the exception.** Partial 23P axons show an orientation rule (−0.056) that full 23P axons don't (−0.002); the difference is −0.054 (SE 0.028, z = −1.9).
  - The calibrated simulation predicts ~0, so if real, this is not a completeness effect. It could reflect which 23P neurons got fully proofread.
  - It is one of 16 uncorrected tests. Treat it as a lead to check, not a finding.

## Revised summary for the preprint

1. Incompleteness *can* bias an inferred wiring rule, but how much depends on the pattern of loss and on each cell type's spatial wiring layout.
   - Cutting axons back uniformly from the far end inflates the RF like-to-like rule by up to 32%. The effect comes from 4P axons, whose RF-similar partners lie near the soma.
   - Random loss of the same amount of cable has no effect.
2. Simulated with MICrONS's actual partial-proofreading profiles, the bias is ~1% (CI within about ±10%).
3. Real partial axons agree with full ones within error, as the calibrated simulation predicts.
4. **So the published like-to-like rule is robust to the reconstruction incompleteness present in MICrONS.**
   - The method (calibrate a cut model on partial vs. full axons, then simulate) gives an error budget for any wiring rule, which mirrors the Fly project.
   - Interareal (V1→HVA) projections are the exception: they are lost first under cutting back from the far end.

## Next

- **Final run on the pinned environment** (`scripts/run_session7_finish.sh`, ~5 h):
  - the main bootstrap with rotated reps
  - calibrated fits, bootstrap and shuffle
  - the within-cell-type radial bootstrap
  - partial vs. full with bootstrap and shuffle
- **Ld (Tweedie) bootstrap:** optional.
- **Sensitivity:** block size; per-cell-type calibration uncertainty (resample the Session 5 profile).
- **23P orientation lead:** do fully vs. partially proofread 23P differ in depth or position, or in tuning?
- **Write the preprint around the revised summary above.**
