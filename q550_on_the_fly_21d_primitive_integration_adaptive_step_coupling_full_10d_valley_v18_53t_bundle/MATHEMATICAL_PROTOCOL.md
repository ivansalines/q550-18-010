# Mathematical protocol — Q550 v18.53t

Let `z(ξ)∈R^11` satisfy the original authenticated Q550 equation `z'=f(z,ξ;κ₀)`, with `z₁₁` the original material clock coordinate. For `Y=z[:10]`, let `q(Y)∈R⁵` be the original embedding, and `F_U(Y,ξ)`, `F_C(Y,ξ)` be the original literal radial-potential and phase-coupling forces. Define the antisymmetric matrix

`B(Y,ξ) = q⊗(F_U+F_C) - (F_U+F_C)⊗q`.

The historical 10D coordinate representation `biv10(B)` uses `sqrt(2)*B_ij` on the ten independent upper-triangular index pairs. The xi-clock source is

`g(z,ξ;κ₀) = (dτ/dξ) biv10(B)`, with `dτ/dξ = ω/sigmoid(z[8])`.

The 21D extension is

```
z' = f(z,ξ;κ₀)
C' = g(z,ξ;κ₀),   C(ξ_initial)=0.
```

**The 11D RHS contains no C dependency.** Thus differences between the physical 11D trajectories of independent 11D and augmented 21D solves are exclusively numerical (including the numerical error norm affecting step acceptance), not a new feedback force.

## Numerical interventions

- **A (11D-only):** DOP853, rtol 3e-12, atol 3e-14, max_step_xi 0.004375; accepted solver steps and 4x-grid dense output saved.
- **B (21D augmented):** DOP853, same physical rtol/atol/max_step; ten additional components use the same rtol and explicit accumulator atol 3e-14. Accepted step mesh may differ from A.
- **C (locked 10D):** On A's frozen numerical dense solution `z_A(ξ)`, solve `C'_locked=g(z_A(ξ),ξ)`, with the same scalar tolerances. No numerical feedback on z_A is possible.
- **D (sampled-source):** Reconstruct `C` from source samples using source PCHIP integration and composite trapezoid as done in v18.53s.

For every variant, `N(a,w)=C(a+3w)-C(a)` is measured in Euclidean 10D and optimized over the same `w` interval as v18.53s. Report minimum norm, `u=N/||N||`, oriented `u∧du` (45 coordinates), directional speed and direct geodesic arclength.

A supplementary check evaluates B's **dense 21D solution C** directly at the optimized endpoints and compares it with cubic Hermite interpolation of grid-sampled C and the source. This separates a second interpolation error from the integrated 21D primitive.

## Diagnostics and interpretation

- `||z_B-z_A||`: physical-state difference introduced by an augmented numerical solve.
- `||g_B-g_A||`: 10D source difference on the same xi grid.
- `||C_B-C_locked||`: accumulator difference between integrated 21D and locked 10D histories; it includes both physical-trajectory and numerical integration differences.
- `||C_locked-C_sourcePCHIP||`: difference between an adaptive integral along one frozen dense 11D trajectory and a source-sampled PCHIP integral.
- `||C_B-C_sourcePCHIP||`: the total numerical-method contrast, not a causal decomposition.
- Accepted step counts, max-step-bound fraction and symmetric nearest-mesh gap quantify step-choice changes.

The stationary radial condition at an interior 10D minimum is `J^T N=0`, **not** `N=0`. A rotating unit direction at such a positive minimum is compatible with ordinary vector geometry and does not certify a physical spin or holonomy.

## Restrictions

All analyses are conditional on the formal clock continuation past `xi=25.6`. Numerical agreement among solver variants is evidence of finite-tolerance robustness, never a rigorous enclosure of the continuum dynamics or of an exact zero. Output descriptions and code comments remain in English.
