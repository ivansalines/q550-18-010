# v18.53m — Mathematical protocol

## 1. Fixed coordinate-free ambient embedding

For a skew 5×5 torque tensor B, the Frobenius-isometric 10D representation is

`vec10(B) = sqrt(2) * (B_ij)_{0<=i<j<=4}`.

For a genuine Q550 solution and material clock τ(ξ), define the **numerically accumulated** torque primitive `C(ξ)=∫_{ξ_min}^ξ vec10(B_U+B_C) dτ`. The implementation reproduces the historical v18.53l trapezoidal primitive and linearly interpolates its accumulated values between native ξ nodes.

The three consecutive residuals are `R_j(o,w)=C(o+jw)-C(o+(j-1)w)` for j=1,2,3. The sum telescopes **exactly** to

`N(o,w)=C(o+3w)-C(o)`.

The lack of internal boundary dependence is a mathematical identity, not a new material law.

## 2. Observer-fixed rank-two plane and two-component root

Compute a rank-2 10D basis `E=[e1,e2]` from the three-tensor SVD at a single predeclared anchor `(o_a,12.8)` and freeze E for the full 2D scan. The governing residual coordinates are

`F(o,w)=E.T N(o,w)`.

A simultaneous 2-component projected zero satisfies `F1=F2=0`. This is distinct from the weaker `min_w ||F(o,w)||` of v18.53l and from the much stronger `N=0` in all ten dimensions. Roots generically are **isolated points**, whereas the finite-tolerance `||F||<=epsilon` can appear as a narrow elongated near-null band or locus.

Three different fixed plane anchors `{13.6,14.4,15.2}` are evaluated independently for each seed and κ₀ rule. An internal constant O(2) change of E cannot change its projected-zero location or its 2×2 Jacobian singular values; changing the rank-2 plane can.

## 3. Analytic boundary attribution and conditioning

Write `s=o`, `e=o+3w`, `g(ξ)=dC/dξ` for the derivative of the **piecewise-linear sampled primitive**. Away from ξ grid knots,

`∂F/∂o = E.T[g(e)-g(s)]`,

`∂F/∂w = 3 E.T g(e)`.

Thus `J=[g_e-g_s, 3 g_e]` where `g_e=E.T g(e)` and `g_s=E.T g(s)`. In particular,

`det J = -3 det[g_s, g_e]`.

The determinant is therefore constrained **solely by the endpoint torque directions in the fixed projected plane**; interior partition boundary changes cannot restore an otherwise singular Jacobian. Report both physical-unit Jacobian singular values and their domain-normalized analogs (scaling o,w by full explored parameter spans). An ill-conditioned Jacobian may turn a narrow near-null into a deceptively precise optimizer solution.

The primitive's segment slope is **exact for that piecewise-linear numerical object**. Separately interpolate actual `B(ξ)*τ'(ξ)` from the real Q550 trace and report disagreement with the discrete-primitive slopes. There is no claim of a closed-form derivative certificate for the underlying continuous ODE at interpolation grid knots.

## 4. Root, false-positive, and stability checks

Use a bounded multistart 2D nonlinear least squares search, report all attempts and deduplicated solutions. A numerical projected root requires successful solver convergence, interior search location and `||F||<1e-8`; a conditioned projected root additionally requires minimum Jacobian singular value >=1e-4, raw Jacobian condition number <=1000, and anchor rank-2 singular gap >=0.003. These are **operational thresholds**, not a proof of existence or uniqueness.

For every candidate compute the **complete** 10D `||N||` and its complement outside the selected plane. Evaluate that same `(o,w)` under the other fixed planes. Re-solve from the same seeds after decimated-trapezoid quadrature; for the four tightened Q550 dynamics re-solve using an independently recomputed torque primitive with the original fixed E.

A candidate at or very near an interpolation knot is explicitly marked; midpoint finite differences may disagree with one-sided derivative as the discrete primitive is only C0 with piecewise-constant derivative.

## 5. Excluded interpretations

No automatic inference of a *physical* flip, new independent axis, physical holonomy, conserved charge, material recurrence, or ten-dimensional bivector closure. These tests locate observer-dependent projected constraints in a formally extended model time trace.
