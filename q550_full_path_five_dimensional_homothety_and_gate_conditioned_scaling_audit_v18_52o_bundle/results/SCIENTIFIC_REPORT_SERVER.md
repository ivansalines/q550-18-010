# Q550 v18.52o — Scientific report

**Status:** AUTHENTIC_FULL_RESOLUTION_POSTPROCESS

## Source and methods

- Authentic v18.52l SHA-256: `d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`
- 4 matched pair labels (original/autonomous × negative/positive), each with amplitudes 1e-4 and 3e-4.
- Common material time via cubic Hermite interpolation using physical 5D dQ and dV; equal array index is NOT assumed to mean equal physical time.
- Interpolator (Hermite vs linear) and naive same-index clocks are explicitly measured as controls.
- Amplitude-dependent gate phase shift is measured separately; sampling gates stay fixed on small amplitude.
- Gate roots are registered on the small-amplitude trajectory alone, with nodal slope m=-0.308 and ±17.04° tilt.
- Baselines include midpoint of successive same-view gates, complementary path, and matched ||q_small|| quartiles.
- All time samples within a path are correlated; no naive time-sample hypothesis test.

## Full-trajectory homothety statistics

| Pair | Median norm ratio | Median relative residual | Max relative residual | Median angle (deg) | c_LS |
|---|---:|---:|---:|---:|---:|
| original/negative | 3.000229559 | 0.000395917689 | 0.000825493083 | 0.0217758 | 3.000183742 |
| original/positive | 2.999774230 | 0.000392709539 | 0.0008181543 | 0.0216129 | 2.999818931 |
| autonomous/negative | 3.000228240 | 0.000394906232 | 0.000822832692 | 0.0217196 | 3.000182656 |
| autonomous/positive | 2.999775521 | 0.000391708885 | 0.000815533128 | 0.0215582 | 2.999820004 |

## Registered gate/midphase comparison

- Registered gate events: 32; phase-matched midpoints: 24.
- Median root relative residual: 0.00037055986; median midphase relative residual: 0.0003744796.
- Midpoint minus adjacent-gate residual: median 7.5731166e-05; positive 16 / 24.
- Direction matters: positive midpoint-minus-gate indicates a gate residual advantage. Negative indicates the opposite.
- Inspect gate_windows_vs_full_path.csv, norm_quartile_conditioned_controls.csv and sensitivity before claiming localized scaling.

## Scope and decision

- A 3:1 ratio is imposed at initialization; persistence over time is a response-linearity test, not new spontaneous autosimilarity.
- Reported 5D homothety is physical only if small and large 5D vectors remain close under the common material clock.
- Equal-looking projected lobes, gate root timing, and exact observer identities alone cannot establish a topological invariant.
- Distinguish a gate-local residual minimum from denominator artifacts and a general nonlinear phase dependence.
- Stage is scientific only when sample_stride=1; nonunit stride reports are DEVELOPMENT_SMOKE_NOT_SCIENTIFIC.

All scripts, source pins, configuration and output SHA-256 manifest are bundled with server sendback.
