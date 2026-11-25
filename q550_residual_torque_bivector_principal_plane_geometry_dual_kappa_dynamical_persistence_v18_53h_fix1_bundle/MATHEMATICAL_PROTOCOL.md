# Mathematical and operational protocol (v18.53h)

## 1. Source hierarchy

The full original v18.51s generator and historical external phase are staged from the SHA-256-authenticated server sendback. The native v18.52l `extended_horizon_10d_states.npz` is also SHA-256-authenticated. v18.53g source-authenticated results and output manifest are pinned separately. No numerical result can substitute for missing original source.

The historical variables are `Y=(log r1,v1,theta1,omega1,log r2,v2,theta2,omega2,z,vz)`, auxiliary `x=sigmoid(z)`, and the saved NPZ's eleventh column `tau` is material time, NOT external phase. The phase `Theta_ext(xi)` comes independently from the staged historical context. Its drive sign comes from the authenticated original ramp clock.

The dynamics in the driving coordinate `xi` uses the authentic `rhs_from_original((Y,Theta_ext),kappa,...)[:10]`, with `d tau/d xi = engine.omega/x`. Each of the two historical seed families supplies one initial condition at `xi=12.8`, common to its κ original and autonomous runs. The historical pre-12.8 path is not replayed separately with the two κ prescriptions.

## 2. Geometric residual

Write `Q in R^5` for the authenticated kinetic embedding and `F_U,F_C in R^5` for the source-exact radial-potential and phase-coupling forces. The physical torque bivectors are `B_U=Q∧F_U` and `B_C=Q∧F_C`. The integrated antisymmetric tensors on an interval I are

`M_U(I) = ∫_I B_U d tau`, `M_C(I) = ∫_I B_C d tau`, `R(I)=M_U(I)+M_C(I)`.

Coordinate-free cancellation is `1-||R||_F/(||M_U||_F+||M_C||_F)`. It is a signed *integrated* cancellation, not a local absolute cancellation; it can depend on the integration interval. The full six-force skew tensor is also measured separately.

Any real skew 5x5 tensor has two orthogonal oriented principal planes with nonnegative singular intensities σ1 ≥ σ2 ≥ 0, plus at least one kernel direction. The 2x2 blocks of a real Schur decomposition yield these planes. Their orientation bases have a sign/gauge freedom. Projector principal angles are invariant under basis rotation, but individual plane identification becomes ill-conditioned near σ1 ≈ σ2. The null direction is not an automatically conserved or dynamically fixed axis.

## 3. Numerical protocols and falsification

- Before new integration: reproduce the archived 53g full-interval radial/phase bivector norms from the exact same native sampling and phase source, within 4e-10.
- Replay shared κ from the same initial state and require max componentwise scaled error ≤2e-4 versus the archived v18.52l trajectory. Report exact error rather than rounding.
- Independent original and autonomous κ IVPs: DOP853, rtol 3e-9, atol 3e-11, max xi step 0.035; both repeated with rtol 3e-11 and atol 3e-13 to check refinement ≤2e-6.
- Every variant checks the material RHS reconstructed by the force decomposition against the authentic source; the native original field is what drives integration.
- The exact counterrotation identity `J=r1²omega1-r2²omega2` is checked along each true trajectory; hinge kinetic moment is not thereby a conserved Noether charge.
- Four disjoint xi windows, as well as the full interval, test whether cancellation and principal-plane structure change through time.
- Compare original and autonomous κ only at the same xi/material context/seed and with identical initial coordinates; phase forcing is historical and prescribed, not recomputed from a new x-dependent clock.

## 4. Scientific classification

A full-pass classification is conditional on source integrity, exact native geometry replay, accepted solver convergence, RHS validation and charge drift. It is not a claim of global bifurcation, full historical κ-separated orbit, state-independent invariant compensation, or canonical hinge Noether charge.

## Fix1 clock-rate regression

Verify the archived material-clock coordinate on each reference native trace using
full-resolution differences, not coarse-stride derivatives:

\[\epsilon_{clock}=\max_i
\frac{\left|\frac{\tau_{i+1}-\tau_i}{\xi_{i+1}-\xi_i}
-\frac12\left(\frac{\omega}{x_i}+\frac{\omega}{x_{i+1}}\right)\right|}
{\max(1,\left|\frac12(\omega/x_i+\omega/x_{i+1})\right|)}.\]

The guard rejects nonmonotone clocks, nonfinite states, or `epsilon_clock > 1e-8`.
The production native input has been checked independently to exhibit
`epsilon_clock < 9e-11` on both reference series. A prior missing variable
caused an audit abort before the initial IVP and is not a physical phenomenon.
