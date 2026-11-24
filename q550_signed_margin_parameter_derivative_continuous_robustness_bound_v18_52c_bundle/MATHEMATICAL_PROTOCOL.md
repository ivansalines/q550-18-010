# Q550 v18.52c — Mathematical protocol

## Domain

The authenticated physical parameter grid comprises 13 distinct kappa values from the 26 nonlinear/tangent branches of v18.52b, across two original physical seed families. The qualified time nodes satisfy `xi >= 0.01` (the first qualified node is `0.01015625`) through `xi = 0.20`. The model's parameter domain is **not** extended in this audit.

## The signed observable

For each seed and physical kappa, the model-derived forced tangent responses produce the source-channel angular budget

\[
m(\xi,\kappa)=\Omega_{\mathrm{mixed}}(\xi,\kappa)+\Omega_{\mathrm{self},v_2,\omega_2}(\xi,\kappa).
\]

Both numerator and denominator are computed in the historical physical Cartesian second plane with the original model channel attribution. A parameter derivative of this margin **at fixed xi** must not be confused with differentiating `min_xi m`, since the minimizing xi can change.

## Three epistemic layers

1. **Authenticated physical samples**: 26 genuine nonlinear/tangent trajectories from v18.52b (13 kappa values x two seeds), with direct hash verification of all five historical result artifacts.
2. **A mathematically continuous surrogate**: piecewise cubic PCHIP in kappa, independently for each qualified xi node. Its exact polynomial extrema and derivative maxima are computed per piece. This establishes properties of the *interpolant only*, not of the nonlinear model.
3. **Independent new physical midpoint probes**: one true nonlinear trajectory plus a new physical 10D forced tangent calculation at each of 12 parameter-gap midpoints per seed (24 total); source jets and angular budgets are recomputed from the historical model. These probes test the surrogate, parameter secants and the observed sign, but still form a discrete set.

## Conditional two-dimensional positivity inequality

Suppose the **true** model function m is differentiable throughout the rectangular domain and independently verified uniform derivative bounds are supplied:

\[
|\partial_\kappa m| \le L_\kappa,\qquad
|\partial_\xi m| \le L_\xi,
\]

plus a uniform absolute evaluation uncertainty bound \(\varepsilon\). Let \(m_*=\min_{i,j}m(\xi_j,\kappa_i)\), and let \(h_\kappa\) and \(h_\xi\) be the maximum spacings. Every domain point is at most half the maximum grid spacing from a sampled node in each direction. By the mean-value theorem,

\[
m(\xi,\kappa)\ge m_*-\frac{L_\kappa h_\kappa}{2}-\frac{L_\xi h_\xi}{2}-\varepsilon.
\]

A positive right side is a **valid conditional mathematical implication**. However, the program's illustrative choice `L = safety_factor * max(observed finite secant, interpolant derivative)` is **not independently proved** to bound the true derivative. In the absence of verified enclosures, any positive numeric right side remains **conditional and non-certifying**. The same applies to the numerical allowance: it is user-assumed, not a validated global error bound.

## Failure and falsifiability

- If any new physical midpoint has a negative margin, record `SAMPLED_PHYSICAL_SIGN_LOSS_FOUND`.
- If solver/tangent reconstruction fails, record `PHYSICAL_MIDPOINT_PROBES_PARTIAL_OR_UNRESOLVED`, preserving logs and a sendback.
- If all physical midpoints are positive, record `PHYSICAL_MIDPOINTS_POSITIVE__CONDITIONAL_CONTINUUM_BOUND_ONLY`.
- **Never** use the category `CONTINUOUS_PHYSICAL_ROBUSTNESS_CERTIFIED` without independent interval derivative and numerical enclosures.

The synthetic test includes a smooth hidden negative interior bump between strictly positive sample endpoints. This rules out the invalid inference that positive samples alone establish continuous stability.
