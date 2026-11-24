# Q550 v18.53b — mathematical protocol

## Fixed angular Noether–d'Alembert charge

`J = r1^2 omega1 - r2^2 omega2`. The 8D quotient coordinates are `(log r1,v1,omega1, log r2,v2,z,vz, Delta)`, where `Delta=theta1+theta2-phi0-Theta_external`. Eliminating `omega2` imposes J **exactly**, not by post hoc projection.

## Constant signed external phase drive

In the historical support coordinate ξ: `dTheta/dξ = omega*s`, `dτ/dξ = omega/x`, thus `dTheta/dτ = s*x`. This is a parametric family of constant-control fields for `s∈[-1,1]`, not a claim that the historical cubic protocol has been re-integrated with a constant s. In the relative phase coordinate `Delta'τ=omega1+omega2-s*x`.

## The stationary branch under test

Set `Delta=π (mod 2π)`, `x=1/2`, `v1=v2=vz=0`. For each `(J,kappa,s)`, solve the *three* algebraic equations

```
0 = omega1^2 - U'(r1)/r1 + kappa*x*r2/r1
0 = omega2^2 - U'(r2)/r2 + kappa*x*r1/r2
0 = omega1 + omega2 - s*x
omega2 = (r1^2 omega1 - J)/r2^2
U'(r) = 2r - 2.4r^2 + 0.8r^3 - e_reserve^2/r^3
```

Unknowns are `(logr1,logr2,omega1)` and `r_i>0`. These equations are independently checked against all eight components of the reduced instantaneous RHS of the source generator.

## Local stability and possible boundaries

For each equilibrium, form `D_u F_reduced` by centered finite differences tangent to fixed J, using two step sizes. The source-based 8D Jacobian must agree with the independently reconstructed local Jacobian and satisfy step convergence. Report all eight complex eigenvalues, the spectral abscissa and the smallest singular value of the **separate algebraic 3×3 equilibrium Jacobian**. Spectral sign crossings and near-singular equilibrium conditions are localized only as *candidate* boundaries. Only if `s != 0` can an auxiliary one-phase-turn multiplier be formed from `exp(lambda * 2π/(omega*abs(s)))`.

For the numerical fold refinement, solve `F(a,s)=0` together with `det D_a F=0` and test the left/right nullvectors `w,v`, the parametric transversality `w^T ∂_s F ≠ 0`, and quadratic curvature `w^T D^2_aa F[v,v] ≠ 0` under three finite-difference scales. Also require one near-zero **real** 8D generator eigenvalue while the other seven stay in the left half-plane. These are numerical nondegeneracy checks and require authentic source confirmation; there is no interval-arithmetic or analytical global proof.

No extrapolation of the source-certified negative-drive Floquet multipliers beyond the v18.53a anchor. At `s=0`, no phase turn exists and the multiplier must be null/undefined.

## Authenticity

Source engine is loaded from the v18.51s archive, pinned to SHA-256 `553775d6...`; v18.53a certified report pinned in `inputs_config.json`. At each stationary state compare the authentic staged `rhs_from_original` against the reconstructed eight-coordinate generator, including `F_11` contribution to `Delta'`. In strict mode any mismatch yields a partial/failure status and nonzero exit.
