# Q550 v18.53l — Mathematical Protocol

## Integrated residuals and exact ambient metric

Use each genuine source pair B_U(xi), B_C(xi) in the kinetic 5D embedding. For width w and origin o, define the three 5×5 skew residuals

R_j(o,w) = integral_{xi=o+(j-1)w}^{o+jw} [B_U(xi)+B_C(xi)] d tau(xi), j=1,2,3.

Represent a skew 5×5 tensor by its ten upper-triangular entries multiplied by sqrt(2). This map preserves the full Frobenius inner product. Let M(o,w) be the 3×10 matrix of these rows and N(o,w)=sum_j R_j (in the 10D embedding).

## Projectors: observer dependence versus declared model transport

Let M=U diag(sigma1,sigma2,sigma3) V.T and P_SVD=V[:,0:2] V[:,0:2].T. This is the rank-2 projector onto the *10D bivector-coordinate* best-fit plane. Its conditioning depends on sigma2-sigma3; it is not a physical 2-plane in the original 5D state.

Anchor P_0=P_SVD(o,w0), w0=12.8. Compare:

1. P_fix(w)=P_0.
2. P_mov(w)=P_SVD(o,w).
3. P_model(w)=L(w) P_0 L(w).T.

For the third case, use the 5D kinetic position q(xi) of the **same integrated physical trajectory**, u=q/||q||, u'=d u/dxi, and the declared minimal radial connection Omega=u' u.T-u u'.T, dR/dxi=Omega R. The code discretizes it by consecutive minimal rotations mapping u_i to u_(i+1), without integrating an arbitrary 10D frame. Relative 5D rotation between origin+3w0 and origin+3w is T=R(xi_endpoint)R(xi_anchor).T. The Frobenius-isometric induced action of B -> T B T.T defines the 10×10 matrix L. This model-dependent construction is one gauge choice and not unique.

For each observer, p(w)=P(w)N(w), and c_a,b=p(w_a)·p(w_b)/(|p(w_a)||p(w_b)|). The model-pullback comparison also computes L(w).T p_model(w) before taking the cosine, so that an apparent ambient rotation due solely to frame transport is separately visible.

## Source-versus-observer decomposition

For two consecutive widths a,b, the identity is exact, without taking a numerical derivative:

Delta(PN) = P_bar Delta N + Delta P N_bar,

P_bar=(P_a+P_b)/2, N_bar=(N_a+N_b)/2.

This is NOT in general a decomposition into orthogonal components, nor an energy partition, because P_bar need not be idempotent and the two terms need not be orthogonal. The code checks the vector identity for all observers and all adjacent width pairs.

## Boundary derivative

With fixed origin o, the sum of three interval integrals is the integral from o to o+3w. Thus:

N'(w)=3 B_tau(o+3w) tau'(o+3w),

with B_tau the original material-clock torque field. As a separate numerical check, use polynomial-exact derivative weights on each uneven 3-node and 5-node local stencil: sum_j c_j (w_j-w_i)^k/k! = 1 if k=1 and 0 otherwise (0<=k<n). Endpoint B_tau and tau' are interpolated from the same numerical IVP trace. This is a quadrature and finite-difference check, not a new physical law.

Additionally check for every adjacent widths a<b:

N(b)-N(a) = integral_{o+3a}^{o+3b} B_tau(xi) d tau(xi),

within numerical quadrature tolerance. Any remaining difference of the projected residual can be assigned algebraically to a varying projector, not necessarily to a material memory effect.

## Conditioning and uncertainty

Require abs(sigma2-sigma3)/sigma2 >= 0.003 for interpretations relying on stable rank-2 projectors; rows violating this are retained and marked as unconditioned. Track projector idempotence, induced SO(10) orthogonality, physical radial transport-map accuracy, and the original v18.53k replay. No numerical sign zero crossing is certified by a two-point negative cosine alone.

## Refusal of overclaim

The transported P_model is a declared *kinematic* connection derived from Q,V, rather than a proven unique physical material transport. A stable fixed-frame inversion is an observation in the invariant Frobenius embedding, not a proof of a hidden mechanical flip. Holonomy would require an actual closed material path and a specified transport rule, not a triangle of observation windows.
