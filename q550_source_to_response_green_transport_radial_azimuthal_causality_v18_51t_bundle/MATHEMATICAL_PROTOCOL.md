# Q550 v18.51t — Mathematical protocol

## Model, scope, and exact provenance

The nonlinear model is the original authenticated 10D driven field from the Q550 genealogy. Its **two** fixed seed states are integrated independently along the original and autonomous κ₀ histories, and its *midpoint* physical state from the v18.51r parent is used for a parameter first variation. The full v18.51s server sendback, itself containing the complete v18.51r ancestry, is the physical source; no substituted force model is introduced. The interval is event-relative `xi ∈ [0,0.2]`, with 513 sections and shared physical drive history for each seed comparison. The two historical κ₀ values are `1/6` and `0.16659961041764054`, and all on-trajectory linearizations use their arithmetic midpoint.

## Source and state variational dynamics

At midpoint state `Y_m(xi)` evaluate the historical analytic state Jacobian and a centered finite-difference parameter derivative at fixed full physical state and fixed drive phase:

\[
J(\xi) = \partial_Y f(Y_m(\xi),\xi;\kappa_m),\qquad
b(\xi)=\partial_{\kappa_0}f(Y_m(\xi),\xi;\kappa_m).
\]

Use two parameter steps, `1e-5` and `2e-5`. Split the **measured** 10D source into three additive pieces in the explicit historical chart:

\[
b=b_{v_2}+b_{\omega_2}+b_{\rm rest},\quad
b_{v_2}=e_5 b_5,\quad b_{\omega_2}=e_7 b_7.
\]

Indices `5` and `7` denote physical velocity variables `v2` and `omega2` in zero-based ordering `(logrho1,v1,theta1,omega1,logrho2,v2,theta2,omega2,logitx,vz)`. Do **not** silently discard other source components.

For each channel integrate a separately forced linear response with **zero** initial parameter sensitivity:

\[
S_j'=J S_j+b_j,\quad S_j(0)=0,
\qquad S_{\rm total}'=J S_{\rm total}+b.
\]

Directly require `S_total = S_v2 + S_omega2 + S_rest` for every section.

## Independent Green representation

Integrate the independent fundamental propagator `Phi'=J Phi`, `Phi(0)=I`, and compute, using interaction-picture quadrature,

\[
S_j(\xi)=\Phi(\xi)\int_0^\xi\Phi(u)^{-1}b_j(u)\,du.
\]

The composite Simpson sum on even section indices and composite 3/8 plus first-edge trapezoid on odd section indices are **independent** from the forced-channel DOP853 integrations. Only even section indices are used for the stringent Green-vs-forced comparison, because they have exact composite-Simpson panel closure. The code records the condition number of `Phi` and aborts on severe ill-conditioning; this is a numerical diagnostic, not interval arithmetic. The model-generated `J,b` nodes are sampled on the nested 129-node and 257-node grids, with cubic interpolants and independent forced-IVP solves to test sampling convergence.

## Physical-coordinate attribution

Map each 10D sensitivity through the analytic immersion differential `D_Y Q`, with `Q=(rho1 cos theta1, rho1 sin theta1, rho2 cos theta2, rho2 sin theta2, sqrt(mx) x)`. At each saved midpoint state project the second oscillator's Cartesian response into the *moving* unit radial basis `(cos theta2,sin theta2)` and azimuthal basis `(-sin theta2,cos theta2)`. Record all six signed channel-to-response coefficients: radial and azimuthal contributions from `v2`, `omega2`, and `rest`. Since the basis is evaluated on the reference midpoint orbit, **sum before taking norms**, and separately record any cancellation.

A nonzero `omega2 -> radial` response, or `v2 -> azimuthal` response, is model-mediated linearized cross-channel transfer. Its magnitude does **not** prove these source components can be independently varied by the underlying physical parameter κ₀ or that they are unique causal mechanisms under an arbitrary coordinate transformation.

## Independent cross-checks and limitations

- Reconstruct the total sensitivity independently against the physical original–autonomous κ₀ finite secant stored in v18.51r.
- Compare the total Green/forced response endpoint against the integrated v18.51s server result.
- Require channel additivity in 10D and in the moving second-plane Cartesian basis.
- Report parameter-step, source-sampling-grid, and Green-quadrature errors separately.
- The statement **no full-state return has been observed** is limited to the tested ξ range; the search for physically meaningful return sections is not decided by this audit.
- Attribution is conditional on the declared explicit state variables and reference chart; no globally objective torque, energy channel, intrinsic memory, or material torsion follows without further work.
