# Session 2 results: usable set for the like-to-like baseline

Run 2026-09-28 with `scripts/build_usable_set.py` at v1822, and revised in Session 3 (new multi-match rule, see step 4). The output is `data/usable_set_v1822_s2.parquet` (694 neurons, not committed). The Session 2 version of this file (662 neurons) is kept as the "lowest residual" sensitivity variant. Full per-stage tables, including the sensitivity variants and provenance, are in `notebooks/02_counts.md`, which the script generates.

## Rules applied

1. **Start:** the Session 1 set of 1,059 neurons (proofread with dendrite ≥ clean, and in `coregistration_manual_v4`).
   - **Root ids:** two proofread cells have `valid_id != pt_root_id`. For both, `chunkedgraph.suggest_latest_roots(valid_id, t_v1822)` returns the v1822 root, which confirms the lineage. Everything joins on the v1822 root. The assessed root is kept as `proofread_root_id`. Only **1 of the 2** is in the 1,059 set (nucleus 296613, 23P). The other has no manual co-registration.
   - **Cell types:** we applied `aibs_metamodel_celltypes_v661_corrections` (queried at v1822, all rows created 2024-02-01), joined on nucleus id. It relabels 6 proofread cells. Only **1** of them is in the set: nucleus 330086 changes from 4P to 23P. The other 5 are inhibitory cells with no manual co-registration.
2. **Excitatory:** we keep rows with `classification_system == excitatory_neuron` after the corrections. This drops the 1 MC.
3. **Match rules:** we drop matches with `score < 0`.
4. **Join + multi-match rule (revised in Session 3):** we keep only matches with a row in `digital_twin_properties_bcm_coreg_v4` (joined on `session, scan_idx, unit_id`). Among those, each neuron keeps **one**: the highest `cc_max`, with ties broken by the higher `score`.
   - The result is stored in `chosen_scan` (for example `s6_scan4`). `n_matches` and `n_scans` count the neuron's matches that pass the score rule and have digital-twin data.
   - The rule picks a single unit rather than a scan, because 41 neuron–scan pairs match two units in the same scan (overlapping fields).
   - Unit keys are unique in the digital-twin table, and the root and nucleus ids agree with the co-registration table on every joined row.
   - The Session 2 rule is kept as a sensitivity variant: lowest `residual` among score ≥ 0 matches, then join.
5. **Quality:** `cc_abs > 0.2 & cc_max > 0.4`.

## Counts (main pipeline)

| stage | n | fully_ext | partially_ext | interareal | none |
|---|---:|---:|---:|---:|---:|
| start | 1,059 | 274 | 656 | 117 | 12 |
| excitatory | 1,058 | 274 | 655 | 117 | 12 |
| match rules | 940 | 244 | 581 | 106 | 9 |
| joined DT | 811 | 204 | 495 | 104 | 8 |
| quality | **694** | **193** | **392** | **101** | **8** |

Final set by cell type × axon level:

| cell_type | fully_ext | partially_ext | interareal | none | All |
|---|---:|---:|---:|---:|---:|
| 23P | 75 | 118 | 59 | 3 | 255 |
| 4P | 92 | 183 | 20 | 2 | 297 |
| 5P-ET | 3 | 19 | 6 | 1 | 29 |
| 5P-IT | 22 | 48 | 13 | 0 | 83 |
| 5P-NP | 0 | 1 | 1 | 0 | 2 |
| 6P-CT | 0 | 10 | 1 | 0 | 11 |
| 6P-IT | 1 | 13 | 1 | 2 | 17 |
| **All** | 193 | 392 | 101 | 8 | 694 |

**Fully extended axons: 193 of 274 survive (70%).**

| cell type | start | match | joined | quality |
|---|---:|---:|---:|---:|
| 4P | 127 | 114 | 96 | 92 |
| 23P | 106 | 92 | 80 | 75 |
| 5P-IT | 29 | 27 | 23 | 22 |
| 5P-ET | 10 | 9 | 4 | 3 |
| 6P-IT | 2 | 2 | 1 | 1 |

## Sensitivity (supplement)

| variant | match rules | joined DT | final | fully_ext in final |
|---|---:|---:|---:|---:|
| **main** (score ≥ 0, max cc_max, QC) | 940 | 811 | **694** | **193** |
| lowest residual (Session 2 rule), QC | 940 | 775 | 662 | 186 |
| Ding match (residual < 20, score > −10), max cc_max, QC | 1,054 | 912 | 780 | 216 |
| main, no QC | 940 | 811 | 811 | 204 |
| lowest residual, no QC | 940 | 775 | 775 | 198 |
| Ding match, no QC | 1,054 | 912 | 912 | 230 |

## Notes and caveats

1. **`residual` does exist in `coregistration_manual_v4`.** The table has `residual`, `score` and `field`. Ding's rule was applied to it directly, with no stand-in column. Every manual match for these neurons has residual < 20 µm (max 19.8), so in practice only the `score > −10` part of the rule removes anything.
2. **The digital-twin join is still the biggest loss (−129 neurons).**
   - Two scans in `coregistration_manual_v4` have no digital-twin rows at all: 5/3 and 7/4. Scan 8/7 is also missing, but it has only 3 rows.
   - Other scans are missing about 5–10% of their units.
   - The Session 3 rule rescues the 36 neurons whose lowest-residual match had no digital-twin data but another passing match did (775 → 811 at the join).
3. **Effect of the cc_max rule.**
   - 145 neurons in the final set had more than one passing match with digital-twin data.
   - Of the 660 neurons in both the new and the Session 2 final sets, 65 now use a different unit.
   - **2 neurons from the Session 2 set drop out** (roots 864691135214986368 and 864691135339103718). Their highest-`cc_max` unit fails `cc_abs > 0.2`, even though another passing unit clears both thresholds.
   - Because the rule selects on `cc_max`, which is also a filter variable, `cc_max` in the final set is biased upward for multi-match neurons. Keep this in mind if `cc_max` is used as a covariate.
4. **The neuron with a changed root** fails the quality filter (cc_abs = 0.04). **The relabeled neuron** has a single match with score −1.4. Neither is in the final set, so these corrections change no final counts. Both are logged in `02_counts.md`.
5. **Axon `none` (8) and interareal (101)** are kept in the file with `axon_level` and an `interareal` flag, so they can be analyzed as separate groups. The completeness axis (partial vs. full) should exclude both.
6. The cell-type/axon-level confound from Session 1 remains. 6P-CT has no fully extended axons, and 6P-IT has only 1.

## Output columns

`root_id, nucleus_id, cell_type, cell_type_v661, cell_type_corrected, axon_level, dendrite_level, interareal, session, scan_idx, unit_id, OSI, DSI, gOSI, gDSI, pref_ori, pref_dir, readout_loc_x, readout_loc_y, cc_abs, cc_max, cc_norm, chosen_scan, n_scans, n_matches, field, residual, score, proofread_root_id, root_changed_since_proofreading`
