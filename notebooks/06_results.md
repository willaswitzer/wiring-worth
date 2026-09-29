# Session 6 results: like-to-like baseline with a shuffled-wiring control

Run 2026-09-28 at v1822.

**Pipeline, in run order:**
1. `scripts/prepare_like_to_like.py` → `data/l2l_pairs_v1822.parquet`, `data/l2l_neurons_v1822.parquet`
2. `.r-env/bin/Rscript scripts/like_to_like.R real` (real fits)
3. `scripts/run_like_to_like_perms.sh` (permutations)
4. `scripts/summarize_like_to_like.py`

**Outputs:**
- Full table (all variants): `notebooks/06_like_to_like_table.md`
- Figure: `figures/06_like_to_like.png` (main variant)
- Raw fits and permutations: `data/like_to_like/` (not committed)

**Runtime:** real fits took 2.7 min (18 variant × similarity × model combinations). The permutations took 75 min, plus a 126-min overnight run (4 workers) that took the main Ld models from 48 to 500 shuffles.

## Data

- **Presynaptic:** the 193 usable neurons with fully extended axons (the reference set).
- **Postsynaptic:** all other usable neurons.
- **Pairs:** 133,749. They carry 3,586 synapses, after excluding synapses with `pre_axon_dist_um > 3` (the dendrite-to-dendrite false detections from Session 5). 71,949 pairs have Ld > 0.
- **Ld:** from the Session 5 relabeled skeletons (`pairs_usable_v1822.parquet`, `labels=synapse`).

## Similarity measures

- **RF-centre distance:** the angle at the eye between the two digital-twin readout locations.
  - The readout coordinates span [−1, 1] across the monitor.
  - Monitor: **31.8 × 56.5 cm (height × width), ASUS PB258Q, 1080 × 1920 px, 15 cm from the eye**. Source: MICrONS Consortium 2025 (Nature, doi 10.1038/s41586-025-08790-w), Methods, "Monitor positioning".
  - It is assumed centred on and normal to the eye, as in Ding et al. 2025, Methods.
  - Median 6.9°, range 0–30°.
- **Orientation difference:** |Δpref_ori| folded to 0–90°. Only pairs where both neurons have gOSI > 0.25 are used: 37,664 pairs, with 107 of the 193 presynaptic neurons tuned.
- Slopes are reported per **10°** of dissimilarity on the log (link) scale. **Negative = like-to-like.**

## Projection types and inclusion

Projection types come from `nucleus_functional_area_assignment` (V1 vs HVA = AL, RL, LM). The inclusion rule is Ding's: > 30 synapses, > 5 presynaptic neurons with synapses, and no neuron contributing > 50%.
- **HVA→V1 (feedback) is excluded everywhere:** 25 synapses from 4 neurons, one of which contributes 72%.
- **HVA→HVA is excluded for orientation** (20 synapses). With extended dendrites only, V1→HVA orientation is excluded too (21 synapses).
- Only 14 of the 193 fully extended axons are in an HVA, so the HVA rows are thin.

## Models

The three Ding et al. 2025 GLMMs, fitted with glmmTMB 1.1.14:

| model | response | family | data |
|---|---|---|---|
| Ld | Ld | Tweedie, log link | all pairs |
| Nsyn | Nsyn | Poisson | all pairs |
| Nsyn/Ld | Nsyn with offset log(Ld) | Poisson | pairs with Ld > 0 |

- **Fixed effects:** similarity × projection type, plus postsynaptic dendrite level.
- **Random intercept:** presynaptic neuron within projection type, `(1 | pre:proj)`. This is Ding's u_k(i,j),i term, a presynaptic-neuron random effect by projection type.
- **Slopes:** `emtrends`. Every real fit converged with a positive-definite Hessian.
- **Nakagawa R²** (conditional / marginal):

| similarity | Ld | Nsyn | Nsyn/Ld |
|---|---|---|---|
| RF distance | 0.97 / 0.50 | 0.31 / 0.17 | 0.03 / 0.003 |
| orientation | 0.97 / 0.44 | 0.25 / 0.14 | 0.02 / 0.004 |

## Shuffled control

- **What is shuffled:** postsynaptic functional properties (readout location, pref_ori and gOSI together) are permuted among all 694 usable neurons within cell type × area, and the models are refitted. Presynaptic properties stay real.
- **p vs shuffle:** two-sided, centred on the shuffled mean: (1 + #{|null − mean| ≥ |real − mean|}) / (N + 1).
- **Convergence:** all 3,128 permuted fits converged.

**Number of permutations (target 1,000):**

| variant | Nsyn and Nsyn/Ld models | Ld models (Tweedie) |
|---|---:|---:|
| main | 228 | 500 |
| + soma distance | 96 | 20 |
| extended dendrites only | 180 | 36 |

**Why not 1,000:** each Tweedie Ld fit over 126k pairs takes about 30 s and peaks at 1.7 GB. On this 8 GB machine, with an unrelated 4 GB virtualization process running, only 2 of those fit at once. The Poisson models ran on 6 workers. The main Ld models were extended to 500 in a separate overnight run (`like_to_like.R perm main tweedie 500 240 4`, 126 min) once the virtualization process was gone.

**Consequence:** the smallest attainable p is 1/501 = 0.002 for the main Ld models, 1/229 = 0.004 for the main Poisson models, and 1/21 = 0.048 and 1/37 = 0.027 for the soma-distance and extended-dendrite Ld models. A p at the floor means the real slope lay outside every shuffle.

## Main results (193 fully extended axons)

| model | similarity | projection | real slope [95% CI] | shuffled mean [95% interval] | p vs shuffle | n perms |
|---|---|---|---|---|---:|---:|
| Ld | RF distance | V1→V1 | −0.894 [−0.917, −0.871] | −0.164 [−0.260, −0.072] | 0.002 | 500 |
| Ld | RF distance | V1→HVA (FF) | −0.354 [−0.434, −0.275] | +0.015 [−0.158, +0.205] | 0.002 | 500 |
| Ld | RF distance | HVA→HVA | −0.856 [−0.960, −0.751] | −0.060 [−0.210, +0.086] | 0.002 | 500 |
| Ld | orientation | V1→V1 | −0.018 [−0.023, −0.013] | −0.000 [−0.013, +0.011] | 0.006 | 500 |
| Ld | orientation | V1→HVA (FF) | +0.005 [−0.020, +0.029] | −0.005 [−0.039, +0.029] | 0.545 | 500 |
| Nsyn | RF distance | V1→V1 | −0.855 [−0.964, −0.746] | −0.252 [−0.402, −0.110] | 0.004 | 228 |
| Nsyn | RF distance | V1→HVA (FF) | +0.004 [−0.357, +0.366] | +0.066 [−0.358, +0.445] | 0.751 | 228 |
| Nsyn | RF distance | HVA→HVA | −0.568 [−1.018, −0.118] | −0.120 [−0.647, +0.332] | 0.070 | 228 |
| Nsyn | orientation | V1→V1 | −0.040 [−0.061, −0.020] | +0.002 [−0.027, +0.028] | 0.004 | 228 |
| Nsyn | orientation | V1→HVA (FF) | −0.117 [−0.235, +0.002] | −0.035 [−0.169, +0.115] | 0.275 | 228 |
| Nsyn/Ld | RF distance | V1→V1 | −0.096 [−0.216, +0.024] | −0.094 [−0.202, +0.013] | 0.952 | 228 |
| Nsyn/Ld | RF distance | V1→HVA (FF) | +0.387 [+0.016, +0.759] | +0.012 [−0.343, +0.365] | 0.044 | 228 |
| Nsyn/Ld | RF distance | HVA→HVA | +0.238 [−0.237, +0.712] | −0.089 [−0.577, +0.393] | 0.179 | 228 |
| Nsyn/Ld | orientation | V1→V1 | −0.021 [−0.041, −0.001] | +0.004 [−0.022, +0.030] | 0.061 | 228 |
| Nsyn/Ld | orientation | V1→HVA (FF) | −0.136 [−0.262, −0.010] | −0.032 [−0.164, +0.108] | 0.170 | 228 |

## Interpretation

1. **The shuffled null is not centred on zero for RF distance.** In V1→V1 the shuffled slopes average −0.164 (Ld), −0.252 (Nsyn) and −0.094 (Nsyn/Ld).
   - Shuffling within cell type × area keeps each group's RF distribution. Different cell types and layers were imaged in different scans, whose RF centres sit in different places, and connection rates differ between groups. That creates an RF–connectivity association *between* groups, which this shuffle preserves.
   - This is exactly why the real slope must be judged against the shuffle, not against zero. Where the Wald CI excludes 0 but the shuffle p is large, the real slope isn't beyond this null.
2. **RF location is like-to-like at the axonal scale, but not beyond the shuffle at the synaptic scale (V1→V1).**
   - Ld and Nsyn fall far more steeply with RF distance than any shuffle: −0.89 vs −0.16 for Ld (beyond all 500 shuffles), and −0.86 vs −0.25 for Nsyn (beyond all 228).
   - Synapses per unit Ld (−0.096) is indistinguishable from the shuffle (−0.094, p = 0.95).
   - Neurons with nearby RFs have more axon–dendrite proximity, and hence more synapses. Given the proximity, they are not converted to synapses at a higher rate.
   - The same holds for Ld in V1→HVA and HVA→HVA.
3. **Orientation is like-to-like in V1→V1 at every scale, weakest at the synaptic scale.**
   - Ld: −0.018 per 10°, p = 0.006 (2 of 500 shuffles as extreme).
   - Nsyn: −0.040, p = 0.004.
   - Nsyn/Ld: −0.021, p = 0.061. Adding soma distance gives p = 0.041, and extended dendrites only gives p = 0.083.
   - Per 10° of orientation difference, synapse count falls by about 4% (e^−0.040). Across the full 90° range that is about a 30% fall.
4. **V1→HVA and HVA→HVA are underpowered** (119 and 89 synapses). Their only nominally significant synaptic-scale result is *anti*-like-to-like for RF distance in V1→HVA (Nsyn/Ld +0.39, p = 0.044; +0.56, p = 0.017 with extended dendrites). This is one of 15 main-variant tests, uncorrected, so don't read it as a finding without more HVA data.

## Sensitivity (full tables in `06_like_to_like_table.md`)

**+ soma-to-soma distance (3D, per 100 µm) as a covariate:**
- The RF effects shrink but stay beyond the shuffle: Ld V1→V1 goes from −0.89 to −0.52 (p = 0.048, the minimum for 20 shuffles), and Nsyn from −0.86 to −0.64 (p = 0.010). So soma distance explains part of the RF effect, not all of it.
- The orientation effects are essentially unchanged: Nsyn −0.036 (p = 0.021), Nsyn/Ld −0.026 (p = 0.041).
- Nsyn/Ld for RF in V1→V1 becomes −0.167 against a shuffled mean of −0.084, still not beyond the shuffle (p = 0.19).

**Extended dendrites only** (103,320 pairs; the orientation analysis is V1→V1 only):
- Same pattern. RF: Ld −0.98 (p = 0.027, the minimum for 36 shuffles), Nsyn −0.97 (p = 0.006), Nsyn/Ld −0.10 (p = 0.97).
- Orientation: Ld −0.022 (p = 0.027), Nsyn −0.047 (p = 0.011), Nsyn/Ld −0.022 (p = 0.083).

## Caveats

- **Permutation counts are below 1,000.** The main Ld models have 500 and the main Poisson models 228. The sensitivity Ld models have only 20 and 36, so their p-values sit at the floor. To go further, rerun `like_to_like.R perm <variant> <track> <n> <min> <workers>`; it resumes. Four workers ran well once the virtualization process was gone.
- **p-values are not corrected** for the 15 main-variant tests.
- **Wald CIs understate uncertainty** relative to the permutation null, because the null's spread is set by the structure in the population, not only by sampling noise.
- **The reference set is fully extended axons only.** The Session 5 finding that partially extended 23P axons are radially truncated therefore doesn't bias these results. It matters if partially extended axons are added later.
- **Not yet used:** Ding's main similarity measures (feature-weight similarity and in-silico signal correlation) aren't in the public digital-twin table.
