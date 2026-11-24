# Q550 v18.52b — Mathematical protocol

## Objective

Determine whether the **signed cross-channel turning margin** measured in v18.52a remains positive after *physical* parameter continuation of the original 10D nonlinear Q550 generator. This is a finite numerical scan, not a rigorous continuous-in-parameter certification.

The parameter is the original `kappa0`, not a source-mixing coefficient and not a change of observation frame. Both original and autonomous seed preparations are retained independently.

## Unchanged model and source convention

For each tested parameter value, integrate the original physical nonlinear state equation using its authenticated solver. Evaluate the original structural state Jacobian and the physical parameter source on that new trajectory:

\[
Y'=f(Y,\xi;\kappa_0),\qquad J=\partial_Y f,\qquad b=\partial_{\kappa_0}f.
\]

For the four forced response columns (the two model-supported channels, the residual and their sum):

\[
S_j'=J S_j+b_j,\quad S_j(0)=0,
\qquad b=b_{v_2}+b_{\omega_2}+b_\mathrm{other}.
\]

The Cartesian immersion and moving angle are inherited **unchanged** from the authenticated v18.52a/v18.51t chain. The channel separation is fixed by the original velocity indices 5 and 7, not by a free transformation of source space.

## Signed turning budget

With source-generated 2D vectors `A` and `B`, and `C` the other-source response, let `W=A+B+C`. At observation \(\xi\), define:

\[
\Omega_{\times}=\frac{A\wedge B'+B\wedge A'}{\|W\|^2},
\qquad
\Omega_\mathrm{self}=\frac{A\wedge A'+B\wedge B'}{\|W\|^2}.
\]

The **signed margin** is

\[
m(\xi,\kappa_0)=\Omega_{\times}+\Omega_\mathrm{self}.
\]

It is the margin that represents the mixed contribution minus the negative self contribution when the self contribution is negative. The full turning rate also contains the residual-source terms. A positive `m` is not by itself identical to positive full turning; both are checked. The signed channel area and the sign of the self contribution are checked separately.

Observe only the qualified section domain \(\xi\in[0.01,0.2]\), to avoid poorly conditioned near-zero response vectors.

## Continuation/radius protocol

The historical center and half-width are

\[
\bar\kappa=\tfrac12(\kappa_\mathrm{original}+\kappa_\mathrm{autonomous}),
\quad h=\tfrac12(\kappa_\mathrm{original}-\kappa_\mathrm{autonomous}).
\]

The six v18.52a trajectories already establish the three anchor values \(\bar\kappa-h\), \(\bar\kappa\), \(\bar\kappa+h\) for both seed preparations. Additional values are \(\bar\kappa\pm rh\), by default \(r\in\{2,4,8,16,32\}\). Each new point uses a genuinely new nonlinear physical path, a freshly evaluated Jacobian/source and an independent forced tangent IVP.

A found sign loss is bracketed in the **dimensionless radius** and refined by new physical integrations. No sign loss implies only a **largest sampled positive radius**, not a mathematically proven stable interval. A solver or source evaluation failure is classified separately. Ordinary numerical tolerances and grid refinement are not interval bounds.

## Independent controls

- Authenticated byte-identical parent v18.52a and v18.51s inputs, complete ancestry.
- Linear addition of forced-channel responses, and recomposition of the angular budget.
- Structural parameter-source finite-difference step consistency at fixed state.
- Extra coarse-jet comparison at the outermost tested parameters.
- Signed-margin, full-turning, oriented-area and negative-self-channel criteria separately.
- A saved per-branch ledger and explicit scientific decision independent of numerical audit pass status.

## Physical scope

This is a **single-parameter continuation** through the historical generator and seed preparations. Neither all physically admissible perturbations nor a global material-memory phenomenon are tested. The channel decomposition and the wedge-angle measure are attached to the model's specified state and Cartesian immersion.
