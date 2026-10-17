# Session 12 results: why the orientation rule looked weaker in partial data; postsynaptic rollback

Run 2026-10-16 (11:28 → 13:33) at v1822, pinned environment, `scripts/run_session12.sh`. All fits converged.
Tables: `notebooks/12_missing_synapses.md`, `notebooks/12_results_tables.md`. V1→V1, Nsyn, slope per 10°.

## 1. The weaker orientation rule comes from which axons were incomplete

- **Missing vs. attached synapses** (`scripts/missing_synapses.py`). At Nov 2025 the synapses not yet attached
  were onto more orientation-similar targets (mean orientation difference 33.7° vs 38.1° attached; Feb 2026: 31.6°
  vs 38.5°). Permuting attachment among each axon's own synapses reproduces the difference (Nov 2025: null −5.1°,
  p = 0.63; within 100 µm path bins p = 0.92; within 50 µm target-distance bins p = 0.74). Within an axon, the missing
  synapses were a random sample.
- **Composition.** 4P axons, which carry the orientation rule, were the least complete: 52%, 59% and 67% of their
  V1→V1 synapses attached at Apr 2023, Nov 2025 and Feb 2026, against 64%, 75% and 95% for 23P.
- **Count-matched random loss** (`scripts/histsub_sets.py`, 20 reps). Each complete axon keeps a random subset of its
  synapses of the size it had at that state.

  | | complete | real history | count-matched random loss, mean [2.5%, 97.5%] |
  |---|---|---|---|
  | orientation, Nov 2025 | −0.040 | −0.028 | −0.027 [−0.036, −0.018] |
  | orientation, Feb 2026 | −0.040 | −0.023 | −0.023 [−0.033, −0.016] |
  | RF, Nov 2025 | −0.856 | −0.962 | −0.883 [−0.925, −0.819] |
  | RF, Feb 2026 | −0.856 | −0.943 | −0.888 [−0.925, −0.844] |

  Random loss with the real per-axon counts reproduces the orientation slope exactly but not the RF inflation.
  The RF inflation needs the spatial pattern of loss within axons; the orientation weakening does not.
- **History by cell type.** Within 4P axons, orientation beyond the distance-matched shuffle stays significant in
  every partial state: −0.052 (p = 0.015), −0.048 (p = 0.025), −0.039 (p = 0.050) at Apr 2023, Nov 2025, Feb 2026;
  complete −0.055 (p = 0.005). The RF inflation is 4P-specific in the real history too: 4P shift −0.210
  [−0.353, −0.088] at Nov 2025, 23P +0.009 [−0.119, +0.125].

## 2. Postsynaptic rollback changes nothing material

`scripts/postsyn_history.py`: a synapse counts only if attached to both somata. 97–99% of the l2l synapses were
already attached postsynaptically (3,477, 3,546, 3,568 of 3,586 at the three states). Every shift and gap moves by
≤ 0.01; e.g. RF shift at Nov 2025 −0.102 [−0.195, −0.018] (pre only: −0.106), RF beyond distance −0.146
(−0.148), orientation shift +0.014 (+0.012).

## 3. Multiple comparisons

Six headline tests Holm-corrected in the preprint (Table 2): five survive (adjusted p 0.03–0.04); the orientation
shift at Nov 2025 does not (0.14).

## Correction to Session 8

The "V1→V1 synapses attached (of 3,586)" column in `08_results.md` and the Figure 1A curve count all 3,586 synapses
onto co-registered neurons; 3,353 of them are V1→V1. Labels fixed; numbers unchanged.
