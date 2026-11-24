# Mathematical protocol — Q550 v18.52o

Let `q_a(tau) ∈ R^5` denote the physical relative 5D trajectory from v18.52l, at kick amplitude `a`, with the same branch seed and sign. Each of four comparison pairs uses `a=1e-4`, `b=3e-4`. The original `tau` arrays differ slightly across amplitudes (by up to ~1.63e-5 at the start for an inherited example). **Never equate array-index correspondence to common material-time correspondence.** We restrict to their physical overlap and interpolate the larger-amplitude physical `q,dq/dtau` using componentwise cubic Hermite interpolation onto small-amplitude `tau` points. A linear-interpolation negative control and same-index negative control quantify the effects of these choices.

## Full physical 5D homothety

At common material time define

- `r(tau) = q_b(tau) - 3 q_a(tau)`;
- `e_rel(tau) = ||r||_5 / [3||q_a||_5]` (when `||q_a||` exceeds a declared floor);
- `e_abs(tau) = ||r||_5`;
- `e_global(tau) = ||r||_5/[3 median_tau ||q_a||_5]` (fixed denominator within each pair);
- `s_norm(tau) = ||q_b||_5 / ||q_a||_5`;
- `s_D(tau) = D_b/(9 D_a) = (s_norm/3)^2`;
- `theta_5(tau) = angle(q_a,q_b)` using stable `atan2(perpendicular,parallel)`;
- `c_LS = integral(q_a·q_b dtau)/integral(||q_a||² dtau)`, approximated on the dense uniform material grid by summed values. [The script uses an unweighted discrete LS sum; source cadence is nearly uniform.]

A 3:1 **norm ratio alone is insufficient**: vector-angle and residual tests must pass too. An imposed 3:1 input amplitude is not an independent observation of a 3:1 dynamical law.

## Gate registration independent of the tested residual

Fixed observational nodal slope is `m=-0.308` (the four source point coordinates are unavailable); fixed paper tilt `alpha=±17.04°`. The common sloped hinge has unit axial direction `e=(1,m,0,0,0)/sqrt(1+m²)` and planar normal `n=(-m,1,0,0,0)/sqrt(1+m²)`. The hidden direction `h` is taken from the first half of the historical small-kick `original_negative_1.0e-04` trajectory, consistent with v18.52m. For a given 5D point:

- `v=q·n`, `z=q·h`;
- `w_± = cos(alpha) v ± sin(alpha) z`.

The **small-amplitude** `q_a(tau)` alone specifies all gate zero-crossing times (piecewise-linear root interpolation); we do not select or shift the gates to minimize `e_rel`. The same physical 5D residual is compared at small-kick roots, midpoints between successive same-view roots, ±0.10/0.25/0.50 `tau` neighborhoods, and opposite-sign observer windows (`w_+ w_- < 0`). The large-amplitude crossings are measured **only as a separate timing-drift diagnostic**, not as the sampling anchors for the hypothesis test.

## Controls and identifiability

1. Compare event/phase samples and window vs off-window medians for both `e_rel` and `e_global`: a small norm should not produce a fake relative-error gate minimum.
2. Stratify gate and off-gate samples by quartile of `||q_a||` to mitigate confounding by amplitude, with actual group sizes reported. Empty gate strata have explicit missing values, not invented medians.
3. Change the candidate planar nodal slope among `0`, `-0.308`, `tan(-17.04°)` and `+0.308`; change transverse tilt among `0°,10°,17.04°,25°,35°`; change the hidden direction among PCA, Q3, Q4, Q5. Observer selectivity must not be assumed from tautological `w_+ ± w_-` projection identities.
4. Analyze the source `D,G` identities, 5D time-derivative consistency, and the difference between Hermite and linear interpolation. A residual at interpolation scale is not identifiable as nonlinear physics.
5. Report even/odd kick composition from matched `+/-` traces: `q_even=(q_++q_-)/2`, `q_odd=(q_+-q_-)/2` at common `tau`. This checks response parity; with two input amplitudes it does NOT uniquely identify an arbitrary-order nonlinear law.

## Decision rules

- **Global homothety:** residual is small along the full path relative to the stated tolerance, not only at chosen gates. Source amplitude makes this an approximate linear-response finding.
- **Gate-selective refinement:** a positive and consistent midpoint-minus-adjacent-gates difference in multiple seed/polarity pairs, also visible in globally normalized error, observer sensitivity, and norm-conditioned controls. The population of 24 adjacent intervals is serially correlated: it must not be passed off as 24 independent physical experiments.
- **No privileged gate:** if gate residual is comparable to or worse than phase controls, or if the conclusion changes with reasonable observers and norm controls, classify as general near-linearity / observer-conditioned view, not special 5D section invariance.
- **No physical-topological claim:** an arbitrary open 5D trace plus selected 2D projections cannot alone certify a topological invariant, tube, cone or double helix.
