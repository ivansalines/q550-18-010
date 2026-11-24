# Q550 v18.51r — Mathematical protocol

## Separation of interventions

Fix a seed `s` and integrate the inherited nonlinear 10D RHS `Y'=f(Y,xi;kappa)` independently for `kappa_-`, `kappa_+`, and `kappa_m=(kappa_-+kappa_+)/2`. Keep all initial conditions identical within a seed, and both distinct historical seeds as separate preparations.

The chronological material-time factor `c` is computed from the *original* nonlinear RHS and the kinetic immersion, with `Q'=c V`; its integrated history is saved. The comparison is not between two simultaneously existing material boundaries.

## Extended horizon, three budgets, fixed sampling

Use the prescribed uniform grid `xi_j=j(0.2/512)` including the historical `0.05` endpoint, and DOP853 at three predeclared resolution levels. No parameter scanning, endpoint selection or rescaling is allowed after observing an apparently favorable shadowing result. The historic forcing design `beta:2:2` is unchanged.

Compare each replay's saved Cartesian Q at `xi=0.05` with v18.51q's prior physical state to independently detect an inconsistent continuation/reconstruction.

## Orthogonal decomposition and phase test

For each seed define `Delta Q = Q(kappa_+)-Q(kappa_-)`. Using the lower-κ reference velocity `w=Q'(kappa_-)`,

\[
\ell=\frac{\Delta Q\cdot w}{w\cdot w},\qquad
r_\perp=\Delta Q-\ell w.
\]

These quantities measure a local Euclidean phase-equivalent projection, not physical torsion or a theorem of global orbit shadowing. Save both magnitudes, the angle differences, and the raw state difference.

Independently minimize a **single constant lag** `d` on the predeclared time window `[0.05,0.195]`:

\[
J(d)=\operatorname{mean}_{\xi_j\in[.05,.195]}
 \|Q_+(\xi_j)-Q_-(\xi_j+d)\|^2,
\qquad d\in[-.004,.004].
\]

This can rule out a constant shift as a sufficient account, but cannot prove an orbit conjugacy.

## Projected return tests

The inherited material coordinate `y8` is logit(x). Define `x=1/(1+exp(-y8))`; numerically find sign changes of `dx/dxi`, and crossings `x(xi)-x(0)` away from `xi=0`. When an equal-x event exists, evaluate the 5D configuration and 10D state difference relative to the initial state. **Never label a turning point, equal-x projection, or an event-relative chart exit as a material return unless full posture and dynamical criteria pass.**

Numerical zero detections by interpolation on a fixed grid are not validated root intervals. A blank return ledger means no events found in this search, not absence on all intermediate times.

## Classification and limitations

The audit tests finite-horizon *numerical* transverse discrepancy and a defined return search. Three solver budgets provide resolution comparisons rather than certified forward-error bounds. There is no theorem of hyperbolic shadowing, no proof of hysteresis, no parameter-dependent twin surface physics and no identification of the toroidal radius R from Q's 5D state immersion.
