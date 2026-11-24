# Q550 v18.52n — Server scientific report

**Four-Point Nodal-Axis Reconstruction & Paired-Hinge Observer Audit**

- Stage: `AUTHENTIC_FULL_POSTPROCESS`; no new ODE integration.
- Authentic v18.52l 5D trace SHA-256: `d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`.
- Nodal user hypothesis: m = -0.308000000, theta_slope = -17.118832 deg.
- Separate observer tilt: +/-17.04000 deg; tan(-tilt)=-0.306494231.
- Angle magnitude difference: 0.078832 deg.
- Axis source: `USER_PRESPECIFIED_SLOPE_FIXED_INTERCEPT`; line intercept 0.00000e+00.
- Four-point coordinates status: `NO_USER_FOUR_POINT_COORDINATES_SUPPLIED`.
- Historical Q2=0 roots: 4; they are a negative control and not the user four.
- Ancestor v18.52m PCA-axis agreement: 0 deg.
- Hidden PCA train/heldout angle: 10.035320 deg.

## Crucial one-trajectory observer comparison

- Parent v18.52m primary Q1-axis gate duration fraction: 0.0204010.
- Reproduced original Q1-axis gate fraction in current frame: 0.0204033 (piecewise-linear exact gate intervals; v18.52m used sampled trapezoidal integration).
- Proposed sloped nodal axis gate fraction: **0.0219956**.
- Proposed nodal hinge average gate fraction over eight trajectories: 0.0219482; mean gate windows: 4.000.

| Label | Opposite-side fraction | Windows | + crossings | - crossings | Cone R² (diagnostic) |
|---|---:|---:|---:|---:|---:|
| original_negative_1.0e-04 | 0.0219956 | 4 | 4 | 4 | 0.0050 |
| original_positive_1.0e-04 | 0.0219257 | 4 | 4 | 4 | 0.0050 |
| original_negative_3.0e-04 | 0.0220658 | 4 | 4 | 4 | 0.0050 |
| original_positive_3.0e-04 | 0.0218742 | 4 | 4 | 4 | 0.0049 |
| autonomous_negative_1.0e-04 | 0.0219387 | 4 | 4 | 4 | 0.0050 |
| autonomous_positive_1.0e-04 | 0.0218898 | 4 | 4 | 4 | 0.0050 |
| autonomous_negative_3.0e-04 | 0.0220087 | 4 | 4 | 4 | 0.0050 |
| autonomous_positive_3.0e-04 | 0.0218866 | 4 | 4 | 4 | 0.0049 |

## Scientific decision

There is **no reconstructed four-point user line** without actual four coordinate pairs; the user-supplied slope remains a prior hypothesis.
The nodal axis and the separate angular hinge tilt are now correctly distinguished and tested against the v18.52m Q1-axis control.
Opposite-side gate fractions are observer diagnostics, not evidence of two physical filaments or invariant section topology.
The intrinsic hinge-perpendicular radius can vary even where the projected curves look paired.

## Data lineage and limits

The entire analysis uses a **single physical 5D trajectory per kick** in matched material time, eight inherited kick labels.
The ancestor v18.52l vector SHA and v18.52m PCA direction were checked. No high-cost nonlinear integration was performed.
Four-point coordinates may be supplied by --four-points; explicit --adopt-four-point-fit replaces the candidate slope and intercept.
Curvature, winding, ellipses, cones and tubes are not promoted to observer-independent dynamical/topological invariants.
