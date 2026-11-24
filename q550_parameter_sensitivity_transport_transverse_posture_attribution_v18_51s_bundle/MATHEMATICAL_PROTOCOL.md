# Q550 v18.51s — Mathematical protocol

## Actual physical system and historical parameter values

The authenticated historical 10D nonautonomous model is

\[Y'=f(Y,\xi;\kappa_0),\quad Y(0)=Y_s(0),\]

where the two initial seed preparations `original` and `autonomous` remain distinct. The two authentic parameters are `kappa_original=1/6` and `kappa_autonomous≈0.166599610418`, and the midpoint is their arithmetic mean. The forcing protocol and event-relative clock are inherited unchanged from v18.51r and its q→p→... genealogy. **No new physical model is introduced**.

## Level A — independently authenticated finite parameter pair

Using the tightly resolved trajectories saved in the server sendback of v18.51r, compute, at all 513 sections to `xi=0.2`,

\[S_Y^{FD}(\xi)=\frac{Y_+(\xi)-Y_-(\xi)}{\kappa_+-\kappa_-},\qquad S_Q^{FD}=\frac{Q_+-Q_-}{\kappa_+-\kappa_-}.\]

The same seed's initial state is held fixed. The derivative of the authentic 5D Cartesian immersion is independently coded as `chart_jacobian`; require `D_yQ(Y_mid)S_Y^FD ≈ S_Q^FD`, quantify nonlinear midpoint curvature and explicitly retain the fact that `FD` is a secant rather than a variational solution.

The physical Cartesian coordinates are

\[Q=(\rho_1\cos\theta_1,\rho_1\sin\theta_1,\rho_2\cos\theta_2,\rho_2\sin\theta_2,\sqrt{m_x}\,x).\]

At each section split the second-plane sensitivity into radial and azimuthal projections. The first-plane and fifth-coordinate shares are saved separately. Compute the local Euclidean phase-equivalent shift and its orthogonal remainder using the velocity of the same *midpoint* trajectory; neither is asserted to be invariant under arbitrary metric changes.

## Level B — independently integrated nonhomogeneous parameter tangent equation

Reauthenticate and stage the original 10D nonlinear generator, including all 19 historical ancestor checks, using the same parent staging procedure as v18.51r. Along its independently integrated midpoint historical trajectory (saved in the authentic r archive), use the original analytic jet:

\[J_Y(\xi)=\partial_Y f(Y_m(\xi),\xi;\kappa_m).\]

Build the parameter source with two fixed finite-difference stencils at the *same state*:

\[b_{h}(\xi)=\frac{f(Y_m(\xi),\xi;\kappa_m+h)-f(Y_m(\xi),\xi;\kappa_m-h)}{2h},\qquad h=2\cdot10^{-5},10^{-5}.\]

Independently check selected columns of the analytic Jacobian against 2-sided state perturbations. Solve the genuine inhomogeneous first-variation equation

\[S_Y'=J_Y S_Y+b_h,\qquad S_Y(0)=0.\]

The state-space jet and the parameter source are sampled on two fixed section grids (129 and 257 nodes), with cubic interpolation between nodes; DOP853 integrates the linearized equation at all 513 saved sections. Compare both interpolation resolutions, the two parameter-source stencils, and the final variational derivative against the *distinct authentic physical κ trajectories* saved in r. `S_Q = (D_yQ)S_Y` is projected on the two physical Cartesian planes and their radial/azimuthal directions. The difference between full `S_Y` and the finite secant exposes nonlinear parameter response, numerical effects and possible chart discrepancies.

## Decision hierarchy

- If the Level A chart identities fail, the archival baseline is invalid and the full test must not be promoted.
- If the model/source differentiation, linear IVP, grid refinement or cross-comparison fails, the full test fails explicitly; the saved Level A baseline remains an independent, correctly labeled secant.
- A positive result is numerical and local to `xi∈[0,0.2]`, the two seed preparations and the historical parameter pair.
- The transverse fraction and plane attributions are Euclidean chart diagnostics, not objective physical stresses, intrinsic torsion, topological winding or globally certified shadowing.
- **Do not** interpret the failure to find an `x` return by v18.51r as proof that no return exists outside the tested horizon.
