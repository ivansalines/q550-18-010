# Q550 v18.53u — Mathematical protocol

## 1. Why this audit must distinguish two different questions

The authenticated model embeds its five physical position coordinates in the **mass-whitened** Cartesian vector

\[
q=(\rho_1\cos\theta_1,\rho_1\sin\theta_1,\rho_2\cos\theta_2,\rho_2\sin\theta_2,\sqrt{m_x}\,x).
\]

Consequently its canonical kinetic energy is \(T=\frac12\dot q^\top\dot q\). **The canonical kinetic metric is \(I_5\) in these coordinates.** The historical bivector ten-vector uses \(b_{ij}=\sqrt2\,B_{ij}\), \(i<j\), so the induced canonical kinetic metric is already \(I_{10}\). The historical Euclidean 10D norm was not an arbitrary metric choice. Finding it again after a coordinate change is an exact covariance identity, not an independent new physical discovery.

## 2. Exact induced metric and coordinate covariance

Given a physically specified positive-definite five-dimensional metric \(G\), the exterior-square metric in the historical ten-component convention is

\[
H_{(ij),(kl)}=G_{ik}G_{jl}-G_{il}G_{jk},\qquad i<j,\;k<l.
\]

For material generalized coordinates \(\eta=(\log\rho_1,\theta_1,\log\rho_2,\theta_2,\operatorname{logit}x)\), the Jacobian \(J=\partial q/\partial\eta\) gives \(G_\eta=J^\top G_qJ\). If \(L=\Lambda^2J\), the same antisymmetric contravariant tensor is \(b_q=L b_\eta\). Thus \(H_\eta=L^\top H_qL\) and \(b_q^\top H_qb_q=b_\eta^\top H_\eta b_\eta\). Numerically verify this on multiple points of every trajectory for the canonical and two noncanonical trial metrics.

The material coordinate kinetic metric for the historical diagonal case is

\[
G_\eta=\operatorname{diag}(\rho_1^2,\rho_1^2,\rho_2^2,\rho_2^2,m_x[x(1-x)]^2).
\]

The dynamic material basis is used only as a **coordinate audit**, not as a license to compare components at different times without transport.

## 3. Second question: genuine metric sensitivity

There is no unique additional positive-definite non-Euclidean kinetic metric forced by the present canonical 5D Lagrangian. To measure possible **model-conditioned** sensitivity without making up a new fundamental metric, construct two frozen kinetic-velocity occupation matrices from authenticated t physical states:

\[
S_E=\frac{\int_E \dot q\dot q^\top\,d\tau}{\operatorname{tr}\left(\int_E\dot q\dot q^\top\,d\tau\right)}.
\]

Intervals: (i) **preformal** \(\xi\le25.6\), (ii) **full_formal** all of the preserved 4x trace. In the original Q550 clock \(d\tau/d\xi=\omega/x\), with an overall \(\omega\) factor cancelling from the normalized covariance. Use the same authentic historical \(\dot q\) and the fixed physical mass \(m_x=1\).

Define a frozen family (\(\lambda=1,3\)):

\[
G_{E,\lambda}=\frac{I_5+\lambda S_E}{\det(I_5+\lambda S_E)^{1/5}},\qquad
H_{E,\lambda}=\Lambda^2G_{E,\lambda}.
\]

The determinant normalization removes an arbitrary global scale. **The parameter \(\lambda\), epoch and volume normalization are declared modeling choices; these are kinetic-occupation-weighted diagnostic metrics, NOT a unique kinetic metric derived as the Hessian of the physical Lagrangian.** They are frozen across the complete valley curve to prevent metric variation and connection terms from being silently conflated with source transport.

## 4. Minima and angular transport

From the authenticated t on-the-fly 21D primitive \(C(\xi)\), build the **historically consistent** cubic Hermite interpolant with nodes \(C_i\) and authentic derivatives \(g_i\). The residual is

\[
N(a,w)=C(a+3w)-C(a).
\]

For each frozen \(H\) solve \(\min_{a,w}\frac12 N^\top HN\) over the identical historical center window and compare metric-dependent minima. Along a locally optimized width curve \(w(a)\), define \(\tilde u=H^{1/2}N/\sqrt{N^\top HN}\), with the **principal symmetric** positive-definite square root fixing a reproducible orthonormal identification. Accumulate direct geodesic angles and the 45-component oriented polygon \(\sum (\theta_i/\sin\theta_i)\tilde u_i\wedge\tilde u_{i+1}\). Distinguish scalar arclength (metric dependent) from componentwise shape (also identification dependent). The canonical choice is \(G=I_5\).

## 5. Inferences and non-inferences

- Canonical kinetic metric = historical Euclidean metric: mathematical identity in the exact mass-whitened embedding.
- Exact coordinate covariance under log-polar pullback: mandatory check, **not** new evidence of physical persistence.
- Frozen kinetic-occupation weights: genuinely different positive-definite metrics conditioned on the model's measured kinetic motion; metric-dependence is interpreted as sensitivity, NOT as a newly discovered force or canonical kinetic law.
- Metric-positive minima are not physical state returns. Bivector residual is not a full 11D distance.
- The far endpoints at \(\xi\simeq52.7\) are in formal clock continuation past \(25.6\).
- No rigorous continuum error bound, unique global minimum, spin, holonomy, or material closure is claimed.
