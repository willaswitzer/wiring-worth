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

## Next

- **Preprint figure:** plot the paired shift slope(f) − slope(1) with bootstrap CIs (radial vs. both nulls), plus the real-minus-shuffle gap. The current figure plots marginal slopes, so its bands overlap more than the paired contrasts warrant.
- **Ld (Tweedie) bootstrap:** 100 samples, overnight.
- **Per-cell-type slopes:** to test the re-weighting account (point 4).
- **The plan's third curve:** real 23P partial vs. full axons, at their observed cable ratio.
- **Sensitivity:** block size (10 and 50 µm); branch-order truncation.
