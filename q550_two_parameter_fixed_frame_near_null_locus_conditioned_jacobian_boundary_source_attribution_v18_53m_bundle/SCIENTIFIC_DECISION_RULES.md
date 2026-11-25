# v18.53m — Preregistered classification rules

- `FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE`: missing original physical parents, manifest mismatch, native source drift, or runtime exception. Partial CSV rows are not a genuine Q550 result.
- `PREFLIGHT_PASS__NO_NEW_Q550_INTEGRATIONS`: historical parents and synthetic tests authenticated, but NO new physical solutions run.
- `AUTHENTIC_TWO_PARAMETER_NO_RESOLVED_PROJECTED_ZERO__NEAR_NULL_ONLY`: all original six Q550 traces and four tightened repeats computed, but none of the three fixed observer planes has a resolved projected zero satisfying the predefined numerical residual and bound criteria.
- `AUTHENTIC_TWO_PARAMETER_PROJECTED_ROOT_CANDIDATES_ONLY__NO_MATERIAL_CLOSURE`: at least one selected fixed 2D plane has an interior numerical projected root; this class does not imply full-net closure or root robustness.

For each candidate distinguish:

1. `solved_projected_root` — numerical 2D F1=F2 root, strictly within predeclared bounds.
2. `conditioned_projected_root` — also passes the SVD-plane rank-gap and raw 2×2 Jacobian conditioning thresholds.
3. `coarse_root_remains_resolved` — recomputed quadrature with every-other original ξ node remains solvable near the previous location (a *perturbation*, not an independent IVP).
4. `refined_root_remains_resolved` — full stricter ODE integration recalculates torque/clock before constrained reoptimization (only four seed/kappa channels).
5. `at_other_anchor_*_norm` — numerical root re-evaluated under other **fixed planes**, not just another coordinate gauge on the same plane.
6. `full_bivector_norm`, `complement_bivector_norm` — essential falsification of an unqualified full 10D closure claim.

A fixed plane's internal SO(2) gauge must leave its projection norm and Jacobian singular values invariant. Distinct *anchor planes* need not agree. The two-endpoint Jacobian identity has to close at roundoff on the piecewise-linear discrete primitive; a difference between interpolated physical torque and discrete-primitive slopes is measured and never silently reclassified as material memory.

**No universal model law or exact mathematical root is certified from floating-point multistart search alone.** If a root is proposed, additional interval bounds, mesh convergence, and independent physically grounded observer constraints would be needed.
