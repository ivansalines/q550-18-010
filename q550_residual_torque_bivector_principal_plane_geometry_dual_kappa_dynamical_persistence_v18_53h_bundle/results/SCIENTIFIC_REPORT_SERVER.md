# Q550 v18.53h — Residual-Torque Principal Planes & Dual-κ Dynamic Persistence

## Scientific decision
`AUDIT_FAIL_NOT_CERTIFIED`

## Actual source availability
Parent: SHA256_AUTHENTICATED; Native: SHA256_AUTHENTICATED

## Scientific scope
The frozen-state v18.53g provenance is authenticated, but does not certify distinct trajectories.
On full execution, two independent native 10D IVPs per initial seed are run at κ-original and κ-autonomous, plus a shared-κ baseline replay; a stricter replay checks convergence.
The external historical phase is preserved from the authenticated context; the initial state at xi=12.8 is common between κ variants. This is NOT a complete re-integration of pre-xi=12.8 historical ramp.
Integrated skew bivectors are decomposed into two oriented principal planes. Principal-plane tracking near spectral degeneracy must be interpreted only at the combined-subspace level.
Bivector compensation score uses signed integrals and is not time-local cancellation. No globally conserved hinge charge or new bifurcation certificate is asserted.

## Metrics
{
  "diagnostic": "not computed"
}
