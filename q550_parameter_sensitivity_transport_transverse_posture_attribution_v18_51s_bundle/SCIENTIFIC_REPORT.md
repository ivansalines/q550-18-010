# Scientific report — Q550 v18.51s

## Reproducible finding at packaging time

The archived v18.51r high/low κ trajectories have been independently inspected and their source/output SHA-256 manifests verified. With fixed seed and κ separation `6.705624902611507e-5`, the original seed's final finite-difference Cartesian parameter sensitivity has norm `0.046482231900830416`. Its second Cartesian plane `(Q3,Q4)` accounts for approximately `99.3966%` of the **squared** sensitivity norm. In this plane, at `xi=0.2`, the radial component is approximately `0.0397701` and the azimuthal component approximately `0.0237886`; both are nonzero. The autonomous seed is nearly identical but remains a separately prepared historical trajectory.

The positional gap grows almost quadratically on the tested short horizon, compatible with initial equal positions and a parameter-dependent acceleration. This statement alone does not identify the underlying forcing or prove a κ law. Reconstructing the Cartesian derivative using the analytic immersion Jacobian gives a maximum absolute difference of about `6.4e-12` from the finite-difference configuration derivative, across the archived samples.

## What is *not* finished locally

The independently integrated physical sensitivity `S'=J S+partial_kappa f` requires the original heavy authenticated historical generator. The complete server run is designed, scripted, and covered by independent synthetic controls; its empirical comparison is pending actual execution on the user's machine. Do not cite Level B as a completed local calculation.

## What the server report should decide

The server computes a direct parameter-source derivative from the original field at fixed physical state, integrates the associated 10D inhomogeneous variational equation and measures agreement with previously integrated κ trajectories. It also quantifies the sensitivity to the differentiation step and the jet interpolation grid. The response's concentration in the second oscillator can then be attributed to a measured transport source rather than inferred from the secant alone.

The words `shadowing`, `transverse`, `posture`, and `force` retain precisely their diagnostic meaning; no global material-return, hysteretic-memory, torsional invariant or interval-valued numerical certificate is implied.
