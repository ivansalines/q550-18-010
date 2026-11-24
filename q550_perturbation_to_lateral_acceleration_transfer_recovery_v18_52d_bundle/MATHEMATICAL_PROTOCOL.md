# Q550 v18.52d — Mathematical Protocol

## Scientific question

Does an in-model perturbation of the historical Q550 nonlinear trajectory predominantly alter acceleration transverse to the material velocity? If a perturbation is imposed only at the initial section while returning to the same governing parameters, does the later perturbed trajectory exhibit any finite-window return relative to its unperturbed counterpart?

These are deliberately separate questions. A permanent `kappa0` change is **not** a removable perturbation. A 98% lateral acceleration difference would not by itself imply recovery or energy dissipation.

## Two source-grounded experiments

**Experiment A — persistent physical parameter separation.** Reuse six genuinely integrated v18.51r trajectories at original / midpoint / autonomous `kappa0`, with both original and autonomous initial seeds. Compare the original and autonomous `kappa0` at the same *elapsed material clock* `tau`, using independent monotone cubic interpolation. Archived velocity time derivatives supply a lightweight check. In server mode the physically defined acceleration is independently reconstructed from the authenticated historical nonlinear 10D generator:

`f(y, xi; kappa0)`, `H(y) = (Q(y), V(y))`, `c = d tau / d xi`,

`a_tau(y, xi; kappa0) = [D_y V(y)] f(y, xi; kappa0) / c(y, xi)`.

The `v18.52c` result (24 independent midpoint parameter probes) is authenticated as provenance for the robustness discussion, but cannot replace Experiment A's true state and acceleration.

For a baseline tangent `e(tau) = V_ref(tau)/||V_ref(tau)||`, define

`delta_a = a_pert(tau) - a_ref(tau)`

`delta_a_parallel = (delta_a dot e) e`

`delta_a_perp = delta_a - delta_a_parallel`

`lateral_fraction = ||delta_a_perp||/||delta_a||`.

The **signed** scalar parallel component is retained. Do not conflate `||P_perp delta_a||` with `||P_perp(pert) a_pert - P_perp(ref) a_ref||`; their difference is a change-of-tangent contribution which is separately recorded. Comparing at fixed tau slightly changes the external drive phase xi. The phase difference is measured and **not** silently erased.

**Experiment B — finite initial transverse velocity kick.** At the historical midpoint `kappa0`, keep the governing nonlinear vector field and drive event identical, and change **only the initial velocity** by `delta V0 = +/- epsilon ||V0|| e_perp` with `epsilon in {1e-4,3e-4}`. The physical five-dimensional Cartesian velocity kick is mapped back to the five historical velocity coordinates via the exact immersion Jacobian. The initial Cartesian position is held fixed. Evolve each resulting nonlinear trajectory independently with the authenticated field, then align to elapsed material tau.

The kick is an admissible initial-state perturbation of this local ODE, but this alone does not certify preservation of any additional off-model global charge or constitutive constraints. The test describes the specified 10D model and parameter domain only.

Measure the separate `delta Q(tau)`, `delta V(tau)`, `delta a(tau)`, and signed/lateral components. Define a dimensionless diagnostic phase distance using a **declared** length scale `L = ||V0|| tau_end`, `D = sqrt((||delta Q||/L)^2+(||delta V||/||V0||)^2)`; the metric is part of the diagnostic and is not privileged by the model. Report `D_end/D_initial` and `D_end/D_peak`, without treating a finite-window ratio as asymptotic recovery.

## Scientific decision constraints

1. **Lateral-dominant acceleration** is classified only when the physical RHS reconstruction agrees with independently differentiated archived velocities within a stated numerical tolerance and the median post-initial transverse share exceeds 50%.
2. **Finite-window contraction** is reported only as a diagnostic. An asymptotic recovery certificate is always `false` for this finite simulation, even if `D_end < D_peak`.
3. A change in clock must not be hidden: matched material-time comparison and fixed external-phase comparison are distinct observational sections.
4. Changes of the tangent projector and of the acceleration must not be conflated.
5. Numerical integration failures, mismatch with the authentic historic trajectory, or kick-coordinate inconsistencies must be declared separately from scientific sign/strength findings.

## Boundaries of inference

The parameter-pair comparison is over `xi in [0,0.20]` only; this is not a rigorous full Q550 orbit or return-section certificate. We do not assume a restoring force, a Lyapunov function, a dissipative mechanism, or continuous parameter-interval validity. A positive v18.52c signed margin does not imply lateral energy absorption.
