# Q550 v18.52q — Counter-Rotating Exact RHS Conservation Audit

**Analytical identity for archived angular RHS:** d/dt(r1² ω1 - r2² ω2) ≡ 0. This is a symbolic result from the explicit angular equations, NOT inferred from small numerical drift.

## Archived nonlinear RHS probes

- Random off-orbit states: 240; charge-derivative max |residual|: 2.887e-15 (physical clock), 3.281e-14 (support clock).
- Anti-rotation equivariance max |ΔF|: 1.021e-14.
- Deliberately unbalanced angular forcing produces nonzero dC/dt: 0.00347081999 (predicted 0.00347081999).

## Authenticated native Q550 traces
- Series: 20; largest C range: 3.175237850427948e-14.

## Authenticated ORIGINAL structural_jet probe (separate from 2026-Aug archived RHS)
- Status: VERIFIED; tested off-trajectory states: 96; max charge derivative: 5.218048215738236e-15.
- Direct original-vs-archived RHS comparisons: 96; max relative mismatch: 0.8656162640278354.

## Interpretation

- Counter-rotation leaves Delta=theta1+theta2-phi0-Theta exactly unchanged. Both angular momenta receive the same torque. Their difference C is an exact first integral in the documented RHS.
- The generalized Noether–d’Alembert interpretation is justified for the angular sector, including radial/auxiliary forces with zero angular work along X=∂theta1−∂theta2; a purely conservative full 11D Hamiltonian is NOT asserted.
- The charge is real-valued, not a topological winding integer.
- A positive finite-horizon collective radius is not a theorem of global no-collapse: this requires an all-future kinetic-energy upper bound and C≠0.
- The invariant is in native absolute polar coordinates, not in the observer 17.04° tilted plane nor in the difference of two trajectory separations.

## Decision
**ARCHIVED_EXACT_IDENTITY_AND_NATIVE_TRACES__ORIGINAL_ENGINE_NOT_VERIFIED**
