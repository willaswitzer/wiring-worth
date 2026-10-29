# Session 14: cell-type connection preference across the proofreading history (auto-generated, v1822)

V1->V1 candidate pairs of the 193 complete axons. log OR = log odds of connecting to a target of type A rather than B; distance-stratified = Mantel-Haenszel over 50 um soma-distance bins; depth x lateral-stratified = Mantel-Haenszel over 50 um bins of the target's signed depth relative to the presynaptic soma x 50 um lateral bins. Shift = state minus complete [95% paired bootstrap over presynaptic axons, 1000 samples]. Random loss = each axon keeps a random subset of its synapses of its Nov 2025 size (20 reps; mean [2.5%, 97.5%]).

## Fraction of complete-state connections present, by target type

| presynaptic   | target   |   connections (complete) |   kept Apr 2023 |   kept Nov 2025 |   kept Feb 2026 |   median soma distance (um) |
|:--------------|:---------|-------------------------:|----------------:|----------------:|----------------:|----------------------------:|
| 23P           | 23P      |                      365 |            0.76 |            0.78 |            0.95 |                          78 |
| 23P           | 4P       |                      158 |            0.52 |            0.68 |            0.89 |                         132 |
| 23P           | 5P-IT    |                       78 |            0.36 |            0.63 |            0.99 |                         235 |
| 23P           | 5P-ET    |                      165 |            0.61 |            0.75 |            0.98 |                         309 |
| 4P            | 23P      |                      675 |            0.41 |            0.47 |            0.57 |                         134 |
| 4P            | 4P       |                      723 |            0.59 |            0.64 |            0.72 |                          82 |
| 4P            | 5P-IT    |                      227 |            0.65 |            0.73 |            0.78 |                         111 |
| 4P            | 5P-ET    |                      142 |            0.62 |            0.76 |            0.8  |                         193 |

## Target-type preference (log odds ratio)

| presynaptic   | contrast     | log OR                     |   complete | Apr 2023 shift       | Nov 2025 shift       | Feb 2026 shift       | random loss (Nov 2025 counts) shift   |
|:--------------|:-------------|:---------------------------|-----------:|:---------------------|:---------------------|:---------------------|:--------------------------------------|
| 23P           | 5P-ET vs 23P | raw                        |       1.44 | -0.27 [-0.47, -0.08] | -0.07 [-0.21, +0.08] | +0.03 [-0.02, +0.13] | -0.04 [-0.09, +0.04]                  |
| 23P           | 5P-ET vs 23P | distance-stratified        |       4.4  | -0.30 [-0.62, +0.12] | -0.23 [-0.46, +0.07] | +0.03 [-0.07, +0.27] | -0.11 [-0.25, -0.01]                  |
| 23P           | 5P-ET vs 23P | depth x lateral-stratified |       1.83 | -1.40 [+nan, -0.09]  | -0.97 [+nan, +0.00]  | -0.19 [-0.79, +0.00] | -0.75 [-1.41, -0.40]                  |
| 23P           | 4P vs 23P    | raw                        |      -1.25 | -0.38 [-0.63, -0.17] | -0.13 [-0.29, +0.03] | -0.07 [-0.19, +0.01] | -0.12 [-0.16, -0.03]                  |
| 23P           | 4P vs 23P    | distance-stratified        |      -0.47 | -0.25 [-0.46, -0.06] | -0.12 [-0.28, +0.04] | -0.10 [-0.24, +0.01] | -0.16 [-0.24, -0.07]                  |
| 23P           | 4P vs 23P    | depth x lateral-stratified |      -0.31 | -0.35 [-0.65, -0.05] | -0.28 [-0.57, +0.01] | -0.17 [-0.46, +0.00] | -0.29 [-0.40, -0.19]                  |
| 4P            | 5P-ET vs 4P  | raw                        |       0.94 | +0.04 [-0.15, +0.19] | +0.16 [+0.02, +0.28] | +0.09 [-0.04, +0.20] | +0.01 [-0.10, +0.09]                  |
| 4P            | 5P-ET vs 4P  | distance-stratified        |       1.38 | +0.50 [+0.32, +0.69] | +0.49 [+0.31, +0.66] | +0.26 [+0.11, +0.42] | +0.05 [-0.06, +0.14]                  |
| 4P            | 5P-ET vs 4P  | depth x lateral-stratified |       0.14 | -0.00 [-0.62, +0.47] | +0.18 [-0.23, +0.58] | +0.18 [-0.26, +0.59] | -0.06 [-0.32, +0.25]                  |
| 4P            | 5P-IT vs 4P  | raw                        |       0.31 | +0.10 [-0.02, +0.21] | +0.13 [+0.03, +0.22] | +0.08 [-0.00, +0.16] | +0.03 [-0.07, +0.10]                  |
| 4P            | 5P-IT vs 4P  | distance-stratified        |       0.35 | +0.17 [+0.04, +0.29] | +0.18 [+0.06, +0.28] | +0.10 [+0.02, +0.18] | +0.04 [-0.06, +0.12]                  |
| 4P            | 5P-IT vs 4P  | depth x lateral-stratified |       0.11 | +0.05 [-0.14, +0.23] | +0.06 [-0.11, +0.22] | +0.09 [-0.07, +0.26] | +0.04 [-0.07, +0.13]                  |
| 4P            | 23P vs 4P    | raw                        |       0.33 | -0.38 [-0.53, -0.23] | -0.33 [-0.48, -0.19] | -0.25 [-0.38, -0.13] | -0.01 [-0.09, +0.08]                  |
| 4P            | 23P vs 4P    | distance-stratified        |       0.62 | -0.18 [-0.34, -0.03] | -0.19 [-0.35, -0.05] | -0.18 [-0.32, -0.05] | +0.01 [-0.08, +0.09]                  |
| 4P            | 23P vs 4P    | depth x lateral-stratified |       0.62 | -0.20 [-0.40, +0.00] | -0.20 [-0.40, -0.00] | -0.18 [-0.41, +0.05] | -0.02 [-0.16, +0.09]                  |
