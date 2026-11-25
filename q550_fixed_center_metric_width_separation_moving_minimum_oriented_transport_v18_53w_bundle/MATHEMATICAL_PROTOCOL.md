# Mathematical protocol — Q550 v18.53w

Let `C(xi)` be the authenticated 10D cumulative bivector primitive from the 21D on-the-fly replay. For origin `a` and width `w`, define

`N(a,w)=C(a+3w)-C(a)`.

The canonical 5D kinetic metric is `I_5`; its induced exterior-square metric is `H0=I_10`. For each epoch (`preformal`, `full_formal`), the occupancy tensor `S` is computed with the **historic** velocity and physical-clock weighting. The sensitivity metric is

`G_lambda=(I+lambda*S)/det(I+lambda*S)^(1/5)`, `H_lambda=Lambda^2 G_lambda`.

These are not new canonical kinetic metrics.

## Factorial design
At each of 204 nodes (`seed`, `kappa_rule`, `epoch`, `lambda`), let `p0=(a0,w0)` be the canonical minimizer and `p1=(a1,w1)` the conditional metric minimizer recorded by the authenticated parent v18.53v.

For `c,b,h` in `{0,1}`:

- `a_i(c) = a_c + offset_i`, with 141 offsets in `[-0.35,0.35]`.
- `w_i(c,b) = argmin_[w0-0.05,w0+0.05] N(a_i(c),w)^T H_b N(a_i(c),w)`.
- `N_i(c,b)=N(a_i(c),w_i(c,b))`.
- `Q(c,b,h)=sum_i theta_i / sin(theta_i) [u_i wedge u_(i+1)]`, with `u_i = sqrt(H_h) N_i / ||sqrt(H_h)N_i||` and `theta_i=acos(u_i dot u_(i+1))`.

`Q000` is the canonical fixed-center baseline. `Q011` is the **actual old v18.53v polygon**: its sampled widths are reloaded from the authenticated v CSV, rather than numerically refitted. `Q111` is the newly computed moving-center polygon.

All vectors `Q` have 45 real components using the same 45 pair ordering.

For factors `F={center,width,metric}`, take all six permutations `pi` and define:

`D_f = (1/6) sum_pi [ Q(S_pi_before_f union {f}) - Q(S_pi_before_f) ]`.

Then `Q111-Q000 = D_center + D_width + D_metric` **identically**, including interactions. We report each `||D_f||`, `D_f dot (Q111-Q000)/||Q111-Q000||^2`, the closure residual, and separate endpoint angles. We **do not sum angles**.

## Residual-at-minimum attribution
Set `X(p,H)=sqrt(H)N(p)` at `p0` and `p1`. Use four variants `X(p0,H0)`, `X(p1,H0)`, `X(p0,H_lambda)`, `X(p1,H_lambda)` for a two-factor Shapley split of vector displacement and, separately, the change in its norm. Even this exact allocation is an identity of comparison, **not a causal/material attribution**.

## Explicit falsification checks
1. Authenticated parent v, u, t and source helper.
2. 204 distinct parent nodes, archived 45D polygon and 141-point width curve for each.
3. `||Q011-Q_parent|| <= 1e-7` (normally close to floating-point precision because archived widths are reused).
4. 45D three-factor vector closure, 10D two-factor vector closure, and scalar norm closure <= 5e-11.
5. Baseline lambda zero yields zero geometry change, up to roundoff.
6. No implied material holonomy, 11D recurrence, observer-independent geometry, or theorem about the continuum solution.
