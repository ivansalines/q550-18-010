# Q550 v18.52z — Mathematical protocol (English)

## I. Authentic input and geometry

The numerical root is initialized from the SHA-256-authenticated **v18.52y `results/report.json`**; it is not fitted from the v18.52w or v18.52x plots. Physical state is `Y=(l1,v1,theta1,w1,l2,v2,theta2,w2,z,vz)`, with `r_i=exp(l_i)` and conserved `J = exp(2*l1)*w1-exp(2*l2)*w2`. The fixed-charge quotient is `u=(l1,v1,w1,l2,v2,z,vz,alpha)`, `alpha=theta1+theta2 mod 2pi`. At fixed `J` and counterrotation gauge `g`, lift by `w2=(exp(2*l1)*w1-J)/exp(2*l2)`, `theta1=(alpha+g)/2`, `theta2=(alpha-g)/2`. This recovers the full 10D state and guarantees all perturbations retain the exact prescribed charge.

The genuine RHS is the authenticated historical `physical.rhs_from_original` staged from the SHA-pinned v18.51s sendback; the clock `ctx['theta']` is the authenticated **external late-return phase**, not `Y[:,10]` (which is material time tau). The stroboscopic period `T=2*pi/omega` must satisfy `Theta(t+T)-Theta(t)=-2*pi` over THREE periods on the actual late return branch (`cubic_clock_sign=-1`). The fixed native kappa is `0.1666331385421536` shared by both seed families. The parent and 10D native input archives are SHA-256 verified; missing or mismatched inputs fail closed.

## II. Iterated relative fixed-point refinement

Let `P_J(u)` be the physical 10D flow by one source-phase period, reduced to the eight-dimensional quotient. Seed Newton from the **actual v18.52y candidate** for each family. The residual is `R(u)=wrap(P_J(u)-u)`, with only the angular sum wrapped. Compute a central 8D Jacobian on the exact `J` leaf, solve `(DP-I) δu=-R` using least squares, enforce a scaled trust cap and backtracking, and record every accepted update or failure. The accepted root must satisfy **an independently recomputed absolute eight-coordinate norm** `||R||<1e-8`, not merely an optimizer's internal or relative error.

## III. Independent 1/2/3-turn closure, source semigroup and gauge

For the candidate `u*`, independently integrate genuine native 10D trajectories for T, 2T and 3T. Verify `||P^n(u*)-u*||` for `n=1,2,3` with thresholds `(1e-8,3e-8,5e-8)`, and compare `P(P(u*))` versus direct 2T, `P(P(P(u*)))` and `P(P^2(u*))` versus direct 3T. Repeat these gates using **an independently tighter DOP853 solve** and check the trajectory charge drift and agreement. The source is periodically driven only on the audited late-return branch; full historical forward–return periodicity is not inferred.

## IV. Conditional Floquet and numerical resolution

ONLY after all closure gates pass, evaluate `DP_J(u*)` by independent central differences in all eight fixed-J directions with steps `(1e-4,3e-5,1e-5)` and a tighter-source solver check. Eigenvalues of `DP_J(u*)` then provide *conditional Floquet multipliers* on the eight-dimensional `J` leaf modulo counterrotation and under the specified periodic external drive. A spectrum is certified only when the normalized matrix differences for both step comparisons and the independent solver are below `0.005`. Stable/attracting within that leaf is asserted only when every multiplier has modulus at least `0.001` below one. If root or spectral convergence fails, retain the raw diagnostics without claiming Floquet certification.

## V. Energy work–dissipation ledger

ONLY upon closure, reintegrate the authentic 10D candidate over one period with dense output and evaluate its energy in the external phase `Theta(xi)` using the authenticated seven model parameters. On the late branch the material-time phase rate is `dTheta/dtau=-x`, with `d tau/d xi=omega/x`. Audit on two quadrature grids

`Delta E = W_drive + W_x - D_rho - D_x`.

A closed relative orbit has `Delta E ≈ 0` and work–dissipation balance over its period, without requiring the instantaneous energy to be conserved. Compare both quadrature resolutions and retain their residuals. This neither proves global Hamiltonian dynamics nor a uniformly positive radial bound as `t → infinity`.

## VI. Scientific boundaries

* A **late-return phase-stroboscopic relative cycle** is not automatically a complete historical forward-and-return cycle.
* A charge conserved by the RHS is not full energy conservation; within-leaf attraction is compatible with dissipation and drive.
* An eigenvalue at an unclosed seed point is NOT a Floquet multiplier.
* If the root is not closed within tolerance, report `NO_RELATIVE_CYCLE_CERTIFICATE` and do not produce spectral or cycle-energy claims.
* Chart norms and singular values are coordinate-dependent; only the fixed-point map eigenvalues are locally coordinate-conjugacy invariant.
