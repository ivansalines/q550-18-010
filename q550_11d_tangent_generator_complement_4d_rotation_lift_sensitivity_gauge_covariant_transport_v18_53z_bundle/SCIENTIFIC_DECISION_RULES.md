# Q550 v18.53z — predefined decision rules

1. **Reject any inference** if source SHA, parent y manifest, or nested x/w/v/u/t chain fails.
2. Recovered historical external phase must produce a max Frobenius-isometric 10D source fitting error under `1e-6` and unit-circle defect under `1e-6`. Do not replace the inferred phase with an invented sinusoid.
3. Jacobian central differences at two stencils must agree within `5e-5` relative Frobenius norm. The passive tau column is exactly zero.
4. The 5D lifts must obey `Jq L=I`. The fixed-Cartesian-velocity lift additionally must satisfy `Dv L=0` to `2e-11`.
5. Any finite-lag 5D tangent map with nonpositive polar determinant or tiny minimum singular value must be rejected; do not silently force an orientation.
6. The complement spin must be skew, orthogonal to the physical velocity tangent and yield an SO(5) frame. The minimum-turn y baseline must replay within `5e-6` degrees before any new attribution.
7. Gauge pullback and induced generator covariance must hold within `5e-8`. Gauge covariance is an identity check, NOT a certificate of material memory.
8. Report BOTH horizontal lifts and the half-horizon local-spin comparison. No complement rotation is model-unique unless a distinct principle selects the lift and the horizon dependence is controlled independently.
9. All results beyond xi=25.6 are formal clock continuation. Positive residual nor oriented transport proves 11D material return or material holonomy.
10. Treat every derived angular difference as a diagnostic of the declared frame and metric, NOT a physical absolute torsion angle.
