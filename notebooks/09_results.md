# Session 9 results: sensitivity of the distance-matched shuffle

Run 2026-10-13, 21:42–22:56, at v1822 on the pinned environment. Script: `scripts/run_distance_sensitivity.sh`. 200 shuffles per condition; every fit converged, with no NA slopes. V1→V1, Nsyn (Poisson), slope per 10°. The p values are two-sided and centred on the shuffled mean (floor 0.005).

## 1. Bin width (pooled presynaptic cell types)

Each cell is real minus distance-matched shuffle (p).

| state | real | 25 µm bins | 50 µm bins (main) | 100 µm bins |
|---|---|---|---|---|
| RF, complete axons | −0.856 | −0.121 (0.035) | −0.140 (0.005) | −0.172 (0.005) |
| RF, history Nov 2025 | −0.962 | −0.141 (0.025) | −0.148 (0.025) | −0.178 (0.010) |
| RF, radial f = 0.3 | −1.132 | −0.083 (0.28) | −0.122 (0.085) | −0.184 (0.035) |
| orientation, complete axons | −0.040 | −0.036 (0.005) | −0.037 (0.005) | −0.037 (0.005) |
| orientation, history Nov 2025 | −0.028 | −0.023 (0.10) | −0.023 (0.10) | −0.026 (0.045) |
| orientation, radial f = 0.3 | −0.057 | −0.036 (0.040) | −0.037 (0.035) | −0.046 (0.015) |

- **RF:** coarser bins absorb less retinotopy, so the gap grows from −0.12 to −0.17. At every bin width, the gap on the Nov 2025 state equals the complete-axon gap within 0.02. The conclusion that real incompleteness leaves the RF preference beyond distance unchanged holds at all three widths.
- **RF, radial truncation:** the f = 0.3 gap varies more with bin width (−0.08 to −0.18). At 25 µm it is smaller than on complete axons, and it is not significant at 25 or 50 µm.
- **Orientation:** the complete-axon gap does not depend on bin width. The plateau-era state is borderline: p = 0.10 at 25 and 50 µm, 0.045 at 100 µm. "Not distinguishable from the distance-matched shuffle" must become "borderline (p = 0.045–0.10 depending on bin width)".

## 2. Within presynaptic cell type (50 µm bins)

| | real, complete | real, f = 0.3 | raw shift | shuffle shift | change in real − shuffle gap | gap, complete (p) | gap, f = 0.3 (p) |
|---|---|---|---|---|---|---|---|
| 4P, RF | −0.783 | −1.113 | −0.330 | −0.323 | −0.007 | −0.141 (0.020) | −0.149 (0.11) |
| 23P, RF | −1.091 | −1.126 | −0.035 | −0.108 | +0.074 | −0.195 (0.050) | −0.121 (0.31) |
| 4P, orientation | −0.058 | −0.076 | −0.018 | −0.016 | −0.002 | −0.055 (0.005) | −0.057 (0.025) |
| 23P, orientation | −0.001 | −0.018 | −0.017 | −0.016 | −0.002 | +0.011 (0.67) | +0.009 (0.77) |

- **The 4P radial inflation is entirely distance.** The distance-matched shuffle shifts by −0.323 under truncation, against −0.330 for the real slope. The RF preference beyond distance in 4P is unchanged (−0.141 to −0.149). This confirms the shared-thalamic-input explanation within 4P (C. Schneider-Mizell).
- **The orientation rule is a 4P rule.** 4P axons show orientation like-to-like beyond distance (−0.055, p = 0.005); 23P axons show none (+0.011, p = 0.67). This matches the Session 7 partial-vs-full table (23P orientation slope −0.00 on full axons).
- **23P RF:** the preference beyond distance is −0.195 (p = 0.050) on complete axons. Its change under truncation (+0.07) is within the shuffle's spread.

## What changes in the preprint

- RF gap beyond distance: "−0.14" becomes "−0.12 to −0.17 depending on bin width; equal on complete and Nov 2025 states at every width".
- Plateau-era orientation: "not distinguishable (p = 0.10)" becomes "borderline (p = 0.045–0.10)".
- 4P explanation: "fits, not tested within cell type" becomes "confirmed within 4P: the shuffle moves with the real slope (−0.32 vs −0.33)".
- New result: the orientation rule is carried by 4P axons; 23P axons show none.

## Not run

- Ld (Tweedie) models under the distance-matched shuffle.
- Distance-matched shuffle at Feb 2026.
- Per-cell-type history states.
