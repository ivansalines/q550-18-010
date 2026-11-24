# Q550 v18.51y — Mathematical Protocol

**Title:** Model-Constrained Source Subspaces & Physical Attribution Admissibility Audit.

## Scope

This audit studies the already authenticated Q550 v18.51t/v18.51w/v18.51x source and response tensors for the two historical seed families on `xi in [0,0.2]`. No nonlinear trajectories or new variational IVPs are solved. It tests mathematical admissibility classes under specific constraints; it does **not** assume independent controllability of the two source components of the *single* kappa derivative.

## 1. Anchored source supports

In the historical 10D state chart, the fixed-parameter source `b = partial_kappa f(Y, kappa, xi)` has an observed component support in four velocity slots. The two labelled channels are

`b_v2 = P_5 b,  b_omega2 = P_7 b`, with `P_i = e_i e_i^T` (zero-based indices). Define `P_other = I-P_5-P_7`. In this fixed physical chart and with these *fixed labels*, the decomposition is uniquely specified by disjoint support and exact reconstruction `b = P_5 b + P_7 b + P_other b`.

For the X mixing `M_eta = [[1-eta,eta],[eta,1-eta]]`, the sum remains unchanged. Unless `eta=0` (and both source components are nonzero), the first mixed source has some component in support `P_7` and the second has some component in support `P_5`; hence **neither obeys the strict historically anchored support criterion**. This is a conditional characterization, not a universal prohibition of all other physically meaningful channel families.

For each source, the support-leak fraction is measured by the norm of the forbidden-support projection divided by the norm of the mixed source. This uses Euclidean coordinates and is explicitly a source-chart diagnostic, not an invariant physical energy fraction.

## 2. Dynamic response is not support isolated

Evaluate at the 257 original model-Jacobian samples

`l_i(xi) = ||(I-P_i) J(xi) e_i||_2 / ||J(xi) e_i||_2`, `i=5,7`.

A nonzero value measures immediate Jacobian coupling into other state components, **not** a transmitted Green energy fraction or a proof of irreversible mixing. The physical source can have strict coordinate support while the transported response spans multiple coordinates.

## 3. Coordinate chart covariance vs source relabelling

Let `Z = T(xi) Y`, where `T` is an invertible, differentiable 10D chart map, here chosen as a pure rotation of the two velocity-coordinate slots by `theta(xi) = 15 deg * sin(pi*xi/0.2)` with fixed original physical channel labels. This is a *coordinate description* of the same dynamic system, not a new force law. Correct tangent quantities obey

`J_Z = T J_Y T^{-1} + T' T^{-1}`,

`b_Z = T b_Y`, `S_Z = T S_Y`, and `P_i^Z = T P_i^Y T^{-1}`.

For each archived forced response channel `S_i`, test

`J_Z S_Z,i + (T b_i) = T (J_Y S_Y,i + b_i) + T' S_Y,i`.

The original physical response is recovered as `S_Y = T^{-1} S_Z`. The chart-invariant observables are understood as `L_Z = L_Y T^{-1}` and therefore the original physical second-plane vector and mutual angle are unchanged. This is a coordinate-covariance identity verified with archived Q550 matrices, not an independent nonlinear replay.

**Negative controls:** (a) projecting rotated sources onto frozen numerical indices instead of transporting `P_i` creates apparent cross-support leakage; (b) omitting `T'T^{-1}` in `J_Z` creates a false tangent-generator defect. The audit measures both.

## 4. Admissibility outcomes and limits

- **Strict fixed-chart support + physical labels:** anchored original decomposition is unique, and nonzero X mixing is excluded from that class.
- **Consistent state chart changes:** can mix numerical slots freely as a representation, provided all sources, projectors, tangent generators and observables transform covariantly. These do not alter the labelled physical decomposition.
- **Arbitrary sum-preserving source rebookkeeping:** mathematically allowed for linear reconstruction; cannot by itself be called an independently admissible physical intervention.
- **Independently actuated channels / constitutive action / material memory:** not established by these records.

The near-right angle from W is meaningful for the original physically labelled source split in the original physical output plane; it is invariant under *passive chart changes applied consistently*, but not under arbitrary source-channel redefinitions as shown by X.
