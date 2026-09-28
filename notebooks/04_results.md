# Session 4 results: skeletons, axon path length, and proximity (Ld)

Run 2026-09-28 at v1822 on the Session 3 usable set: 694 neurons, of which 585 are presynaptic (fully or partially extended axons). Scripts, in run order:
1. `scripts/strategy_counts.py`
2. `scripts/fetch_skeletons.py`
3. `scripts/compute_proximity.py`, which generates `notebooks/04_proximity_counts.md`

Outputs (not committed):
- `data/skeletons/` — one `.npz` per root, plus the fetch and proximity logs
- `data/synapses_usable_v1822.parquet` — the synapse table with 3 new columns
- `data/pairs_usable_v1822.parquet` — the pair table

**Runtime:** the skeleton fetch took 2 min 58 s (694 roots, 6 threads, batches of 60, 0 failures). The proximity computation took 1 min 55 s. The full pipeline runs in about 5 min.

## 1. Proofreading strategy columns (usable set)

`proofreading_status_and_strategy` has two strategy columns, `strategy_axon` and `strategy_dendrite`, plus two boolean columns, `status_axon` and `status_dendrite`. The usable set's `axon_level` *is* `strategy_axon`, so that crosstab is diagonal. The informative table is dendrite strategy × axon level:

| strategy_dendrite | fully_ext | partially_ext | interareal | none | All |
|---|---:|---:|---:|---:|---:|
| dendrite_clean | 2 | 73 | 48 | 0 | 123 |
| dendrite_extended | 191 | 319 | 53 | 8 | 571 |

`status_axon` is True for all neurons except the 8 with axon `none`. `status_dendrite` is True for all 694.

**How were the partially extended axons extended?** The table itself does not record this. The table's description points to the strategy manifest (microns-explorer.org/manifests/mm3-proofreading), which says:
- **partially_extended:** "The axon was extended outward from the soma, following each branch to its termination." There is no length threshold and no target-directed rule.
- **fully_extended:** the same pass, and then "every endpoint was identified, manually inspected, and extended again if possible."

So partial and full extension differ by a second pass that rechecked endpoints, not by which branches were followed. The manifest lists `axon_column_truncated` (branches cut at column boundaries), but no neuron in this version has it.

The median resampled axon length per neuron is shorter for partially extended axons:
- 4P: 6.2 mm partial vs. 7.7 mm full
- 23P: 5.4 vs. 7.3
- 5P-IT: 7.4 vs. 9.0

Separately, dendrite completeness varies with axon level: 73 of the 392 partially extended neurons have only `dendrite_clean` dendrites, against 2 of 193 fully extended ones. As postsynaptic cells, `dendrite_clean` neurons may have truncated dendrite tips, which lowers their Ld. `post_dendrite_level` is in the pairs table so this can be controlled for.

## 2. Skeletons

- **Source:** the CAVE skeleton service, `skeleton_version=4`. These are pcg_skel skeletons from the level-2 cache, with SWC compartment labels (1 soma, 2 axon, 3 dendrite) and vertices in nm, rooted at the soma. Root ids are immutable, so they are the v1822 objects. pcg_skel was not needed.
- **Validation:** the root vertex equals the soma point. Postsynaptic sites land on dendrite vertices.
- **Result:** 694 of 694 fetched with **0 failures**. The median dendrite length is 3.42 mm.
- **Two presynaptic skeletons have no axon vertices**: 864691137200012353 and 864691135886393456, both 4P with fully extended axons. The whole skeleton is labeled dendrite, so the service's axon/dendrite split failed. Their 54 synapses have no path length, and their 1,386 pairs have `Ld_um = NaN`. This is the only loss.

## 3. Axon path length of synapses

- **Method:** each non-autapse synapse is mapped to the nearest resampled point on its presynaptic axon. Its path length from the soma along the skeleton is saved as `pre_axon_path_um`.
- **Distribution:** median 286 µm (10th–90th percentile: 155–584 µm, max 2.1 mm).
- `pre_axon_dist_um` stores the distance from the synapse to that nearest axon point, for quality filtering. The median is 0.44 µm.
- **5.3% of synapses (414) are more than 3 µm from the labeled axon.** For 297 of them, the nearest skeleton vertex (median 1.2 µm away) is labeled *dendrite*, so an axon piece was mislabeled. Their path lengths are unreliable. Filter with `pre_axon_dist_um <= 3`.
- Autapses have NaN in all three new columns.

## 4. Proximity (Ld), following Ding et al. 2025

- **Method:**
  - Each edge is cut into ≤ 1 µm pieces, represented by the piece's midpoint and its length.
  - Presynaptic axon pieces and postsynaptic dendrite pieces within 5 µm are found with `cKDTree.query_ball_tree`.
  - Ld is the summed length of the dendrite pieces in proximity.
  - A synapse is in a proximity if it is within 3 µm of any proximity point (axon or dendrite) of its own pre × post pair.
- **Pairs:** 405,405 (585 × 693). Of these, 404,019 have Ld, and **172,689 have Ld > 0** (median Ld 31 µm among them).
- **Synapses outside any proximity: 334 of 7,805 (4.3%).** 303 of those 334 are among the synapses that sit more than 3 µm from the labeled axon. So the misses mostly come from compartment-labeling errors, not from the proximity method.
- 49 connected pairs have Ld = 0, for the same reason.

## 5. Pairs table and synapse density

`data/pairs_usable_v1822.parquet` has these columns: `pre_root_id, post_root_id, Ld_um, Nsyn, Nsyn_in_prox, pre_cell_type, pre_axon_level, post_cell_type, post_axon_level, post_dendrite_level`. Nsyn counts all non-autapse synapses; Nsyn_in_prox counts only those within a proximity.

Synapses per mm of Ld, by pre cell type. The main density is Nsyn_in_prox / ΣLd:

| pre cell type | pre neurons | ΣLd (mm) | Nsyn in prox | syn / mm Ld |
|---|---:|---:|---:|---:|
| 23P | 193 | 3,404 | 2,452 | 0.72 |
| 4P | 273 | 3,623 | 4,022 | 1.11 |
| 5P-IT | 70 | 868 | 776 | 0.89 |
| 5P-ET | 22 | 114 | 80 | 0.70 |
| 6P-IT | 14 | 77 | 101 | 1.32 |
| 6P-CT | 10 | 40 | 33 | 0.82 |
| 5P-NP | 1 | 1.6 | 7 | 4.44 |
| **All** | 583 | 8,128 | 7,471 | 0.92 |

- The 4P pre neuron count excludes the 2 neurons with no axon vertices.
- 5P-NP is a single neuron with only 1.6 mm of Ld, so its density isn't meaningful.
- Within each cell type, density barely depends on axon level. For 23P it is 0.74 (full) vs. 0.70 (partial); for 4P, 1.08 vs. 1.15; for 5P-IT, 0.89 vs. 0.90.
- So partially extended axons have fewer connections mainly because they have less axon length and therefore less Ld, not because they form fewer synapses per unit of proximity.
- The full breakdown is in `04_proximity_counts.md`.

## Caveats

1. **Skeleton compartment labels are the main error source.** They cause the 2 neurons with no axon, the 5.3% of synapses far from the labeled axon, and most of the 4.3% of synapses outside any proximity. Mislabeled axon pieces also count as that neuron's "dendrite" when it is postsynaptic, which inflates its Ld slightly. One fix would be to relabel axon vs. dendrite using presynaptic synapse locations, for example with pcg_skel's synapse-based axon split.
2. **The shuffled-wiring control (CLAUDE.md) is not included yet.** It isn't informative for these per-cell-type densities: shuffling synapses among a presynaptic neuron's pairs keeps its total synapse count and total Ld, so the aggregate density is unchanged. The control becomes essential once like-to-like is computed from `Nsyn` and `Ld_um`. At that point it should be reported alongside the result.
3. The postsynaptic population is only the usable set, so Ld and densities are within that subset.
