# Session 1 results: access + proofreading inventory

Run 2026-09-22 with `notebooks/01_setup_and_explore.py`. Parquet outputs are in `data/*_v1822.parquet` (not committed).

## Pinned version

**1822**, materialized 2026-06-27 14:05 UTC. This was the newest of the 7 public versions (117, 943, 1300, 1507, 1621, 1718, 1822). It is recorded in `config.py`.

## Tables used

| Role | Table |
|---|---|
| Proofreading level | `proofreading_status_and_strategy` (2,334 rows) |
| Cell type | `aibs_metamodel_celltypes_v661` (93,839 rows with a root id; 90,467 after dropping multi-soma roots) |
| Co-registration (primary) | `coregistration_manual_v4` |
| Co-registration (automated, compared) | `coregistration_auto_phase3_fwd_v2`, `apl_functional_coreg_vess_fwd`, `coregistration_auto_phase3_fwd_apl_vess_combined_v2` |
| Functional properties (listed only) | `digital_twin_properties_bcm_coreg_v4`, `digital_twin_properties_bcm_coreg_auto_phase3_fwd_v2`, `digital_twin_properties_bcm_coreg_apl_vess_fwd` |
| Functional area (listed only) | `nucleus_functional_area_assignment` |

## Proofreading status (all 2,334 cells)

| strategy_axon \ strategy_dendrite | dendrite_clean | dendrite_extended | none | All |
|---|---:|---:|---:|---:|
| axon_fully_extended | 9 | 618 | 35 | 662 |
| axon_interareal | 59 | 61 | 0 | 120 |
| axon_partially_extended | 256 | 1,217 | 38 | 1,511 |
| none | 3 | 38 | 0 | 41 |
| **All** | 327 | 1,934 | 73 | 2,334 |

- Rows with `valid_id != pt_root_id` (root changed since assessment): **2**
- Duplicate root ids: 0

## Cell type × axon level (all proofread cells, v661 labels)

| cell_type | fully_ext | interareal | partially_ext | none | All |
|---|---:|---:|---:|---:|---:|
| 23P | 169 | 72 | 295 | 3 | 539 |
| 4P | 209 | 21 | 343 | 2 | 575 |
| 5P-ET | 15 | 7 | 43 | 5 | 70 |
| 5P-IT | 60 | 14 | 121 | 2 | 197 |
| 5P-NP | 0 | 1 | 11 | 2 | 14 |
| 6P-CT | 0 | 2 | 159 | 5 | 166 |
| 6P-IT | 4 | 3 | 211 | 17 | 235 |
| BC | 103 | 0 | 108 | 2 | 213 |
| BPC | 2 | 0 | 83 | 0 | 85 |
| MC | 63 | 0 | 49 | 1 | 113 |
| NGC | 2 | 0 | 22 | 0 | 24 |
| (no type) | 35 | 0 | 66 | 2 | 103 |
| **All** | 662 | 120 | 1,511 | 41 | 2,334 |

## Proofread (dendrite ≥ clean) + co-registered

There are 2,261 proofread cells in total. Counts below are unique neurons.

| Co-registration table | proofread + coreg | axon partial or full | axon full | axon partial |
|---|---:|---:|---:|---:|
| **coregistration_manual_v4** | **1,059** | **930** | **274** | **656** |
| coregistration_auto_phase3_fwd_v2 | 1,359 | 1,224 | 401 | 823 |
| apl_functional_coreg_vess_fwd | 1,452 | 1,314 | 432 | 882 |
| coregistration_auto_phase3_fwd_apl_vess_combined_v2 | 1,437 | 1,301 | 430 | 871 |

Using the manual table, the remaining 129 neurons are 117 axon_interareal and 12 with axon none.

**Cell type × axon level, proofread + manual co-registered (n = 1,059):**

| cell_type | fully_ext | interareal | partially_ext | none | All |
|---|---:|---:|---:|---:|---:|
| 23P | 106 | 72 | 199 | 3 | 380 |
| 4P | 127 | 21 | 300 | 2 | 450 |
| 5P-ET | 10 | 7 | 30 | 1 | 48 |
| 5P-IT | 29 | 14 | 74 | 0 | 117 |
| 5P-NP | 0 | 1 | 1 | 0 | 2 |
| 6P-CT | 0 | 1 | 22 | 3 | 26 |
| 6P-IT | 2 | 1 | 29 | 3 | 35 |
| MC | 0 | 0 | 1 | 0 | 1 |
| **All** | 274 | 117 | 656 | 12 | 1,059 |

### How the co-registration tables differ

This comparison is restricted to proofread neurons.

| Automated table | in both | auto only | manual only | shared neurons with an identical functional unit |
|---|---:|---:|---:|---:|
| auto_phase3_fwd_v2 | 1,036 | 323 | 23 | 952 / 1,036 (92%) |
| apl_functional_coreg_vess_fwd | 1,020 | 432 | 39 | 932 / 1,020 (91%) |
| combined_v2 | 1,039 | 398 | 20 | 970 / 1,039 (93%) |

- **Manual (`_v4`)** matches were made by hand. It has the fewest neurons (15,433 unique) and is the highest-confidence table.
- **auto_phase3_fwd_v2** uses the Allen phase-3 spline with linear assignment. By design it matches **excitatory neurons only**.
- **apl_functional_coreg_vess_fwd** uses a B-spline registration on blood vessels (JHU APL) and includes inhibitory cells.
- **combined_v2** merges auto_phase3 and APL. Where they disagree, it keeps the pair with the lower residual.
- For about 8% of the neurons that both a manual and an automated table include, the two assign different functional units. The automated matches need a residual/score filter before we use them.

## Functional-property tables (not yet analyzed)

The three `digital_twin_properties_bcm_coreg_*` tables share the same schema. Each is tied to one co-registration table: v4 = manual, auto_phase3_fwd_v2, or apl_vess_fwd.

- Columns: `session, scan_idx, unit_id, pref_ori, pref_dir, gOSI, gDSI, OSI, DSI, cc_abs, cc_max, cc_norm, readout_loc_x, readout_loc_y`, plus nucleus/root id and position columns.
- Size: `_v4` has 15,780 rows (13,087 roots), `_auto_phase3_fwd_v2` has 40,562 (24,889), and `_apl_vess_fwd` has 40,746 (25,822).
- The table owners ask for citations of Wang et al. 2025 (10.1038/s41586-025-08829-y), Ding, Fahey, Papadopoulos et al. 2025 (10.1038/s41586-025-08840-3), and MICrONS Consortium 2025 (10.1038/s41586-025-08790-w).
- `nucleus_functional_area_assignment` has 144,120 rows. Columns: `tag` (V1/AL/RL/LM) and `value` (distance to the area boundary in µm).

## Surprises and caveats

1. **Almost no inhibitory neurons are co-registered in the manual table.** There are 435 proofread inhibitory cells (BC/BPC/MC/NGC), but only **1** is in `coregistration_manual_v4` and 0 are in auto_phase3. The APL vessel table has **140** and combined_v2 has 86. The like-to-like analysis with manual co-registration is therefore **excitatory only**. Any E→I or I→E functional analysis would depend on the APL automated matches.
2. **The axon levels are confounded with cell type.** 6P-CT (0 full / 159 partial) and 6P-IT (4 / 211) are almost never fully extended. Among excitatory cells, full extension is concentrated in 23P and 4P. MC and BC are the most often fully extended overall. This is the confound to control for in the December analysis.
3. **`axon_column_truncated` is not present in this version.** The only "biased" level is `axon_interareal` (120 cells, 117 of them co-registered). These should be kept separate from the partial/full completeness axis.
4. **103 proofread cells have no v661 cell type.** 99 are missing from the v661 table and 4 are multi-soma roots. 73 of them have dendrite `none`, which means they are probably axon-only reconstructions with no soma in the volume. They drop out of the proofread (dendrite ≥ clean) set anyway.
5. **Cell-type corrections:** `aibs_metamodel_celltypes_v661_corrections` relabels 6 proofread cells. For now we use the uncorrected v661 labels.
6. **Manual co-registration quality:** among the 1,059 usable neurons, 212 of the match rows have a negative score (a competing nucleus was closer). 288 neurons match more than one functional unit because they were imaged in several scans. We need rules for both before modeling: a score threshold, and a choice between picking the best unit or averaging across units.
7. **2 proofread cells have a root id that changed since assessment.** They are still in the table and still need a decision.
