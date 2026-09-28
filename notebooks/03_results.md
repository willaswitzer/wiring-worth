# Session 3 results: synapses among the usable set

Run 2026-09-28 with `scripts/pull_synapses.py` at v1822, on the revised usable set (694 neurons; see `02_results.md`). The output is `data/synapses_usable_v1822.parquet` (not committed). Full tables are in `notebooks/03_synapse_counts.md`, which the script generates.

## What was pulled

- **Table:** `synapses_pni_2`, the datastack's default synapse table (337M rows at v1822).
- **Pre:** the 585 usable neurons with an `axon_fully_extended` (193) or `axon_partially_extended` (392) axon.
- **Post:** all 694 usable neurons, including the interareal and axon-`none` groups.
- The query was filtered on both root columns at v1822, in chunks of 25 presynaptic roots. Each chunk was checked for server-side truncation.
- **Columns:** `synapse_id, pre_root_id, post_root_id, size, ctr_x_nm, ctr_y_nm, ctr_z_nm`.
  - `size` is the table's `size` column, in voxels as stored.
  - The position is the synapse center point (`ctr_pt_position`), in nm.
- **Result:** 14,803 synapses, of which 6,944 are autapses and 7,859 are between different neurons.

## Connected pairs (pre cell type × pre axon level)

A connected pair is an ordered (pre, post) pair with at least one synapse, excluding autapses.

| pre cell type | fully_ext | partially_ext | All |
|---|---:|---:|---:|
| 23P | 961 | 1,254 | 2,215 |
| 4P | 1,936 | 1,814 | 3,750 |
| 5P-ET | 27 | 76 | 103 |
| 5P-IT | 320 | 416 | 736 |
| 5P-NP | 0 | 7 | 7 |
| 6P-CT | 0 | 40 | 40 |
| 6P-IT | 19 | 77 | 96 |
| **All** | **3,263** | **3,684** | **6,947** |

- Each connected pair has 1.13 synapses on average (median 1, max 10).
- 567 of the 585 presynaptic neurons contact at least one usable neuron.
- Per-group connection probabilities are in `03_synapse_counts.md`. They are connected pairs divided by n_pre × 693. Overall the probability is 0.017. For 4P, it is 0.030 with a fully extended axon and 0.014 with a partially extended one.

## Caveats

1. **Autapses are detection artifacts and are excluded from all pair counts.** They are 47% of the rows. Their median size is 460 voxels, against 4,964 for synapses between different neurons, and 584 of the 585 presynaptic neurons have some. They are kept in the parquet. Filter with `pre_root_id != post_root_id`.
2. **These are descriptive counts, not a wiring result.** Fully extended axons do show more connections per neuron, but that is expected from axon length alone. Per `CLAUDE.md`, any like-to-like result built on this table must be reported next to the shuffled-wiring control.
3. The postsynaptic population is only the usable set: neurons that pass co-registration and functional QC. Connection probabilities are therefore within this subset, not general cortical rates.
