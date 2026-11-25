# Mathematical protocol — v18.53v

With a frozen trajectory-dependent positive semidefinite velocity occupation matrix \(S\) (trace normalized to one), define

\[
A_\lambda=I+\lambda S,\quad G_\lambda=\det(A_\lambda)^{-1/5}A_\lambda,\quad H_\lambda=\wedge^2 G_\lambda,
\qquad 0\le\lambda\le 4.
\]

The exterior-power metric uses the historical 5D pair order \(i<j\), and
\(H_{ij,kl}=G_{ik}G_{jl}-G_{il}G_{jk}\).

The accumulated bivector residual is

\[
N(a,w)=C(a+3w)-C(a),\qquad E_\lambda(a,w)=\tfrac12N^T H_\lambda N.
\]

Here \(C\) is the authenticated 10D primitive, modeled by a cubic Hermite spline constrained at sampled nodes by the stored model-generated 10D source. The canonical historical model kinetic metric is exactly \(I_5\) in the mass-whitened Cartesian embedding. These noncanonical frozen \(G_\lambda\) are sensitivity probes only.

The analytic Jacobian is

\[
J=[g(b)-g(a),3g(b)],\quad b=a+3w,\ g=C',
\]

and its exact spline-based Hessian is

\[
\nabla^2_{p}E=J^T H J+\begin{pmatrix}
(g'(b)-g'(a))^THN&3g'(b)^THN\\
3g'(b)^THN&9g'(b)^THN
\end{pmatrix},\quad p=(a,w).
\]

Metric derivative uses

\[
G'=(\det A)^{-1/5}\bigl[S-\tfrac15\operatorname{tr}(A^{-1}S)A\bigr],
\]

with \(H'\) obtained by differentiating the induced-metric minors. For an internal stationary nondegenerate minimum, the *implicit slope* is

\[
\frac{dp_*}{d\lambda}=-(\nabla^2_p E)^{-1}J^TH'N,
\qquad
\frac{d\|N\|_{H}}{d\lambda}=\frac{N^TH'N}{2\|N\|_H}.
\]

Compare these with finite secants from neighboring independent lambda nodes; numerical agreement is not a separate physical certificate. Scale the Hessian with \(D=\mathrm{diag}(0.35,0.05)\) and report eigenvalues of \(D\nabla^2E D\), since raw coordinate condition numbers depend on the units of origin and width.

At every node also reoptimize width along a fixed origin scan of 141 samples, reconstruct the 10D unit direction in the \(H^{1/2}\)-whitened frame, and compute the integrated 45D geodesic polygon. Angles between polygons under changing \(\lambda\) use the **principal symmetric SPD square-root identification** and do not constitute physical holonomy.

Screens for internal positive-Hessian minima, alternative multi-start competitors, branch shifts and parent-u regression are prespecified in `inputs_config.json`. A report of regular numerical continuation refers strictly to the tested family and data interpolation.
