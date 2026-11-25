# Predefined decision rules (English)

1. Reject any unauthenticated parent, substituted source, incomplete IVP count, charge tolerance failure, missing ancestor files, or invalid numerical result. A failure is `FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE`.
2. Require all six parent minima to replay in new genuinely integrated cumulative bivectors to within `minimum_replay_norm_abs_tol`.
3. If all replay and finite differentiation tests pass, classify `AUTHENTIC_ORIENTED_10D_VALLEY_AND_BOUNDARY_GEOMETRY__NO_MATERIAL_HOLONOMY` regardless of whether turning rate is appreciable. Report actual values without a hard-coded prediction.
4. If a residual norm falls below 1e-7, label it only `NUMERICAL_CLOSURE_SCREEN__PHYSICAL_RETURN_NOT_CERTIFIED` and explicitly require full physical-state tests. Do not silently relabel an integrated-bivector near-return as full closure.
5. Missing or nondefined normalized directions at zeros are explicitly marked; never divide by tiny r or infer a sign flip from undefined coordinates.
6. Mark finite-difference/primitive interpolation knot crossings, width optimization at bounds, and any dependence on formal xi>25.6. Report adverse sensitivity findings even when earlier checks pass.
7. Report the endpoint-gradient identity for the interpolated primitive separately from comparison against samples of the physical source.
