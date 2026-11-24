# Q550 v18.51w — Mathematical Protocol

**Title:** Channel-Vector Angle, Oriented-Area & Observer-Critical Compensation Audit

## Source and domain

Use the authenticated v18.51v server sendback (including the v18.51u/t/s/r genealogy). There are two seed families (`original`, `autonomous`) and 513 matching sections on the domain `xi in [0, 0.20]`, each with an increasing physical clock `tau(xi)`. This version does not integrate any new nonlinear or variational IVP; it analyzes the saved Green-channel response vectors.

At each section, let `A = response(v2)`, `B = response(omega2)`, `O = response(other)`, and `S = A+B+O` in the **physical Cartesian second plane**, or in a naturally aligned orthonormal radial--azimuthal frame. The channel splitting refers to the chosen model coordinates, and is not itself a generally invariant physical source decomposition.

## Rotational invariants

For nonzero `A` and `B`, compute

- `D(xi) = A · B` (signed dot product),
- `H(xi) = det[A,B] = A_r B_a - A_a B_r` (signed oriented area),
- `Theta(xi) = atan2(H,D)`, an oriented mutual angle in `(-pi, pi]`,
- `d = D/(||A|| ||B||)` and `h = H/(||A|| ||B||)`, which satisfy `d^2+h^2=1`.

Under a passive **proper** orthonormal change of observation frame `R(xi) in SO(2)`, `D`, `H`, `Theta`, channel norms and the vector `S` (after transforming back) remain unchanged. An orientation-reversing reflection changes the sign of `H`, while preserving `D`; do not call signed area an `O(2)` invariant.

An angle cannot be assigned at the initial section where both source responses are zero. The audit explicitly masks this sample, and preserves the undefined value in the arrays rather than silently substituting zero degrees.

## Observer-critical axes

Define an arbitrary passive observer rotation `psi` relative to the **natural radial axis**. The observer's projected radial component of any vector `v=(v_r,v_a)` is

`p_psi(v)=cos(psi)*v_r + sin(psi)*v_a`.

For `||v||>0`, a projected zero occurs on the unoriented observer axis

`psi_crit(v) = wrap_pi(atan2(-v_r,v_a))`, modulo `pi`.

The audit follows the three zero axes `psi_crit(A)`, `psi_crit(B)`, and `psi_crit(S)` at every valid section. In a moving observer rotated by `gamma(xi)`, the predicted critical angle relative to that observer is `psi_crit(v)-gamma`, modulo `pi`. All root projections and the observer covariance identity are checked numerically.

For any observer axis, the displayed cancellation is

`C(psi)=1-|p_psi(A)+p_psi(B)+p_psi(O)| / (|p_psi(A)|+|p_psi(B)|+|p_psi(O)|)`.

This quantity is undefined when the denominator vanishes, and is frame dependent. `C=1` at an axis perpendicular to `S` does **not** imply `S=0`, or a material return / closure event. A root of `p_psi(B)` is likewise a zero projection, not a vanishing channel vector.

## Controls

1. CRC and SHA-256 for the parent executable genealogy, all original V source files and outputs.
2. Independent T second-plane Cartesian response compared with V physical laboratory response.
3. Source channel additivity and norm preservation.
4. Dot/area invariance under V lab coordinates and two genuinely time-dependent, proper observer families.
5. Signed area reversal for a reflected basis (negative control).
6. Exact zero projection at predicted critical axes, with modulo-pi covariance of those axes.
7. Unit-circle identity for normalized dot and area; independent synthetic geometry tests.
8. Explicit handling of the zero-amplitude origin; both seeds studied.
9. Distinguish frozen-history source channels, geometric observers and physical-time relabeling from a constitutive-memory certificate.

## Scientific scope

All inferences are numerical diagnostics within the archived finite interval. Strong geometric invariance here is invariance **under proper passive rotations of the physical second plane at a fixed source-channel decomposition**. It does not by itself prove physical actuator causality, objective material memory, or a return-section closure.
