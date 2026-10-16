# Session 11: truncation at one absolute radius (auto-generated, v1822)

V1->V1, slope per 10°. Shift = slope(level) − slope(complete) [95% paired bootstrap]. Gap = real − shuffled mean (p, number of shuffles).

## Cable and synapses kept

Column axis (x, z) = (728, 851) um, the median of the 193 presynaptic somata. 28 axons have somata more than 100 um from it and are left complete under `cyl` (14 of them in V1).

V1 presynaptic axons:

|                  |   median_cable_kept |   synapses_kept |   synapses_full |   synapses_kept_frac |   median cable kept, 23P |   median cable kept, 4P |
|:-----------------|--------------------:|----------------:|----------------:|---------------------:|-------------------------:|------------------------:|
| ('cyl', 100)     |               0.38  |            2360 |            3472 |                0.68  |                    0.438 |                   0.375 |
| ('cyl', 150)     |               0.588 |            2994 |            3472 |                0.862 |                    0.649 |                   0.584 |
| ('cyl', 200)     |               0.755 |            3220 |            3472 |                0.927 |                    0.779 |                   0.752 |
| ('cyl', 300)     |               0.936 |            3340 |            3472 |                0.962 |                    0.932 |                   0.94  |
| ('pathabs', 150) |               0.105 |             292 |            3472 |                0.084 |                    0.107 |                   0.105 |
| ('pathabs', 250) |               0.309 |            1418 |            3472 |                0.408 |                    0.304 |                   0.32  |
| ('pathabs', 400) |               0.636 |            2607 |            3472 |                0.751 |                    0.587 |                   0.68  |
| ('pathabs', 600) |               0.896 |            3166 |            3472 |                0.912 |                    0.855 |                   0.923 |


## pathabs: same path radius R from soma (um)

| similarity   | model   |   same path radius R from soma (um) |   slope | shift                   | gap, cell-type shuffle      | gap, distance-matched shuffle   |
|:-------------|:--------|------------------------------------:|--------:|:------------------------|:----------------------------|:--------------------------------|
| orientation  | Nsyn    |                                 600 |  -0.041 | -0.001 [-0.006, +0.004] |                             |                                 |
| orientation  | Nsyn    |                                 400 |  -0.05  | -0.009 [-0.022, -0.000] | -0.052 (p = 0.005, n = 200) | -0.043 (p = 0.005, n = 200)     |
| orientation  | Nsyn    |                                 250 |  -0.08  | -0.040 [-0.068, -0.018] | -0.084 (p = 0.005, n = 200) | -0.063 (p = 0.010, n = 200)     |
| orientation  | Nsyn    |                                 150 |  -0.088 | -0.047 [-0.110, +0.018] |                             |                                 |
| orientation  | Nsyn/Ld |                                 600 |  -0.022 | -0.000 [-0.005, +0.003] |                             |                                 |
| orientation  | Nsyn/Ld |                                 400 |  -0.029 | -0.008 [-0.019, +0.001] | -0.033 (p = 0.025, n = 200) | -0.036 (p = 0.005, n = 200)     |
| orientation  | Nsyn/Ld |                                 250 |  -0.056 | -0.035 [-0.060, -0.015] | -0.061 (p = 0.005, n = 200) | -0.056 (p = 0.010, n = 200)     |
| orientation  | Nsyn/Ld |                                 150 |  -0.058 | -0.037 [-0.096, +0.026] |                             |                                 |
| RF distance  | Nsyn    |                                 600 |  -0.89  | -0.034 [-0.067, -0.002] |                             |                                 |
| RF distance  | Nsyn    |                                 400 |  -0.954 | -0.098 [-0.168, -0.038] | -0.687 (p = 0.005, n = 200) | -0.111 (p = 0.035, n = 200)     |
| RF distance  | Nsyn    |                                 250 |  -1.107 | -0.251 [-0.376, -0.117] | -0.763 (p = 0.005, n = 200) | -0.049 (p = 0.502, n = 200)     |
| RF distance  | Nsyn    |                                 150 |  -1.192 | -0.337 [-0.644, +0.035] |                             |                                 |
| RF distance  | Nsyn/Ld |                                 600 |  -0.115 | -0.019 [-0.060, +0.023] |                             |                                 |
| RF distance  | Nsyn/Ld |                                 400 |  -0.131 | -0.035 [-0.121, +0.027] | -0.029 (p = 0.642, n = 200) | -0.069 (p = 0.274, n = 200)     |
| RF distance  | Nsyn/Ld |                                 250 |  -0.173 | -0.077 [-0.232, +0.080] | -0.034 (p = 0.657, n = 200) | +0.012 (p = 0.886, n = 200)     |
| RF distance  | Nsyn/Ld |                                 150 |  -0.228 | -0.133 [-0.476, +0.294] |                             |                                 |

## cyl: column cylinder radius Rc (um)

| similarity   | model   |   column cylinder radius Rc (um) |   slope | shift                   | gap, cell-type shuffle      | gap, distance-matched shuffle   |
|:-------------|:--------|---------------------------------:|--------:|:------------------------|:----------------------------|:--------------------------------|
| orientation  | Nsyn    |                              300 |  -0.04  | +0.000 [-0.001, +0.002] |                             |                                 |
| orientation  | Nsyn    |                              200 |  -0.044 | -0.004 [-0.008, +0.000] | -0.046 (p = 0.005, n = 200) | -0.039 (p = 0.005, n = 200)     |
| orientation  | Nsyn    |                              150 |  -0.049 | -0.009 [-0.017, +0.001] |                             |                                 |
| orientation  | Nsyn    |                              100 |  -0.048 | -0.007 [-0.019, +0.008] | -0.051 (p = 0.005, n = 200) | -0.036 (p = 0.020, n = 200)     |
| orientation  | Nsyn/Ld |                              300 |  -0.021 | +0.000 [-0.001, +0.002] |                             |                                 |
| orientation  | Nsyn/Ld |                              200 |  -0.024 | -0.003 [-0.007, +0.001] | -0.029 (p = 0.035, n = 200) | -0.032 (p = 0.010, n = 200)     |
| orientation  | Nsyn/Ld |                              150 |  -0.028 | -0.007 [-0.015, +0.003] |                             |                                 |
| orientation  | Nsyn/Ld |                              100 |  -0.025 | -0.004 [-0.016, +0.010] | -0.031 (p = 0.055, n = 200) | -0.030 (p = 0.035, n = 200)     |
| RF distance  | Nsyn    |                              300 |  -0.879 | -0.024 [-0.037, -0.012] |                             |                                 |
| RF distance  | Nsyn    |                              200 |  -0.916 | -0.060 [-0.096, -0.026] | -0.648 (p = 0.005, n = 200) | -0.130 (p = 0.005, n = 200)     |
| RF distance  | Nsyn    |                              150 |  -0.959 | -0.103 [-0.156, -0.051] |                             |                                 |
| RF distance  | Nsyn    |                              100 |  -1.03  | -0.174 [-0.272, -0.084] | -0.724 (p = 0.005, n = 200) | -0.149 (p = 0.010, n = 200)     |
| RF distance  | Nsyn/Ld |                              300 |  -0.109 | -0.014 [-0.031, +0.000] |                             |                                 |
| RF distance  | Nsyn/Ld |                              200 |  -0.111 | -0.016 [-0.053, +0.022] | -0.013 (p = 0.781, n = 200) | -0.075 (p = 0.199, n = 200)     |
| RF distance  | Nsyn/Ld |                              150 |  -0.125 | -0.030 [-0.087, +0.029] |                             |                                 |
| RF distance  | Nsyn/Ld |                              100 |  -0.183 | -0.087 [-0.211, +0.011] | -0.061 (p = 0.348, n = 200) | -0.095 (p = 0.080, n = 200)     |

## Complete axons (reference)

| similarity   | model   |   slope | gap, cell-type shuffle      | gap, distance-matched shuffle   |
|:-------------|:--------|--------:|:----------------------------|:--------------------------------|
| RF distance  | Nsyn    |  -0.856 | -0.597 (p = 0.005, n = 200) | -0.140 (p = 0.005, n = 200)     |
| RF distance  | Nsyn/Ld |  -0.096 | -0.004 (p = 0.945, n = 200) | -0.089 (p = 0.124, n = 200)     |
| orientation  | Nsyn    |  -0.04  | -0.043 (p = 0.005, n = 200) | -0.037 (p = 0.005, n = 200)     |
| orientation  | Nsyn/Ld |  -0.021 | -0.026 (p = 0.055, n = 200) | -0.030 (p = 0.020, n = 200)     |
