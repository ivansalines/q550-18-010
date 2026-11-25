# Mathematical Protocol — v18.53p

Let `C : ξ ↦ R^{10}` be the cumulative bivector primitive exported by the authenticated v18.53o server run. Define for window origin `o` and common width `w`:

`N(o,w) = C(o + 3w) - C(o); r = ||N||; u = N/r` where `r>0`.

Along a candidate width curve `w(o)`, for `e=o+3w(o)` and `g=C'` in differentiable regions,

`dN/do = g(e) * [1 + 3w'(o)] - g(o)`

`dr/do = u · (dN/do)`

`du/do = [dN/do - u (dr/do)] / r`.

The chain-rule angular speed is `||du/do||`, equal to `||u ∧ du/do||` in the ambient Euclidean bivector coordinate metric. The instantaneous oriented form has `10 choose 2 = 45` components.

For an ordered sampled curve `u_i` on the unit sphere `S^9`, its **polygonal great-circle length** is

`L_poly = Σ 2 atan2(||u_{i+1}-u_i||, ||u_{i+1}+u_i||)`.

The corresponding oriented 45D segment-form integral is `(θ_i/sin θ_i) (u_i ∧ u_{i+1})`, interpreted continuously at `θ_i=0`. This is *not* a material holonomy and does not require choosing an SVD gauge.

The **piecewise-linear primitive** has discontinuous derivative at trace knots. The exact bounded minimizer of `||N(o,w)||²` is determined by enumerating all trace segments intersecting `e∈[o+3w_min,o+3w_max]` and minimizing a quadratic along each segment. This avoids local optimizer fallacies but also exposes nonsmooth argmin paths and possible ties. The width path's ordinary derivative need not exist at switching points. Its finite differences can bias the chain-rule tangent-integral estimator.

A **PCHIP primitive** provides a smoother interpolation control over precisely the same sampled cumulative data. It must not be represented as an independent dynamical solution or as a correction of physical equations. The local scalar objective is scanned for bracketed minimizers. Report interpolant sensitivity without selecting whichever result seems to best support a physical claim.

Compare three angular lengths: direct `L_poly`, trapezoidal integration of chain-rule speed using sampled `w'`, and trapezoidal integration using PCHIP derivative of the sampled width. Use a fourth length from finite differences of the unit direction. The direct polygonal length converges to the variation of the *sampled* direction path, not automatically to a unique smooth-ODE arclength.

Measure total width variation, kink proximity of both time boundaries, boundary contact, projected source vs cumulative-interpolant derivative, and the scalar radial identity. The previous 0.0098–0.0099 rad difference is the numerical target; no threshold certifies a physical phase transition.
