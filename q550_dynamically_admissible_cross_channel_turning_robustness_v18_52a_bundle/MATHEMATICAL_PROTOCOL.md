# Mathematical protocol — Q550 v18.52a

## Historical physical family

Let the authentic generator be `Y' = f(Y, theta_drive(xi); kappa0)`, with the same original dynamical constraints and drive protocol as the v18.51r nonlinear replay. Test three historical values:

\[
\kappa_a=0.16659961041764054,\quad
\kappa_m=(\kappa_a+1/6)/2,\quad
\kappa_o=1/6.
\]

For each seed and each historical kappa, use its own preintegrated physical trajectory `Y_kappa(xi)` from v18.51r, rather than interpolating one trajectory to another. At every sampled node recompute

\[
J_\kappa(\xi) = \partial_Y f(Y_\kappa(\xi),\kappa),
\qquad b_\kappa(\xi) = \partial_\kappa f(Y_\kappa(\xi),\kappa).
\]

The source derivative is central-differenced in *kappa at fixed physical state* using steps `h=1e-5` and `2h` for a consistency check; the 10x10 Jacobian comes from the structural analytic jet.

## Model-fixed channel transport

The canonical channel split, defined by physical source slots, is

\[
b_{v_2}=e_5e_5^\top b,\quad
b_{\omega_2}=e_7e_7^\top b,\quad
b_{\rm other}=b-b_{v_2}-b_{\omega_2}.
\]

With identical zero initial parameter sensitivity, independently solve each nonhomogeneous tangent equation

\[
S_j'=J_\kappa S_j+b_j,\qquad S_j(0)=0.
\]

Push these vectors through the analytical Cartesian immersion Jacobian `D_Y Q`, take the second Cartesian plane `(Q_3,Q_4)`, and write its vector responses `A`, `B`, `C` for the three channels. Their sum is the full response `R=A+B+C`.

## Signed turning budget

With the oriented 2D wedge `u wedge v = u_x v_y - u_y v_x`, define

\[
\Omega_R=\frac{R\wedge R'}{\|R\|^2}.
\]

The expansion separates self terms `A wedge A'`, `B wedge B'`, cross term `A wedge B'+B wedge A'`, and all terms involving `C`. All four are divided by the *same* `||R||^2`. The cross-dominance margin is

\[
m(\xi,\kappa)=\frac{A\wedge B'+B\wedge A' + A\wedge A'+B\wedge B'}{\|R\|^2}.
\]

Report the minimum signed margin, positive-turning sample count, signed endpoint budget, vector angle and oriented area across every historical physical branch and both seeds, using `xi >= 0.01` to avoid undefined turning near the zero-response origin.

No claim of strict sign between sampled sections is made without an interval bound. No automatic presumption of robustness is built into the test: if `m <= 0` on a branch, it is an empirical counterexample to that specific stronger condition.

## Controls

The numerical audit independently checks source support in the original state coordinates; decomposition/recomposition in the 10D sensitivity; angular wedge closure; variation of the parameter-difference step; half-density vs full-density jet nodes at the midpoint; and baseline agreement with the v18.51z server-produced responses. The optional fundamental-matrix/Green comparison is from the v18.51t reference implementation.

The three-branch audit may support a *conditional* finding over the specific historical kappa0 band. It does not establish a metric-independent, channel-choice-independent, or all-perturbations theorem.
