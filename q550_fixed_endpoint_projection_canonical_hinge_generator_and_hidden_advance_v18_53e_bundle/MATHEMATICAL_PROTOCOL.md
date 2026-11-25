# Q550 v18.53e — Mathematical protocol

## A. Fixed endpoints, passive observers, and hidden-advance identity

Let `q(t)=deltaQ(t)` be a 5D **difference between two independently integrated Q550 trajectories** (kick minus reference) at matched material time, not an absolute orbit. Its time derivative is `v(t)=deltaV(t)`. Choose the *same* two times `t_a<t_b` once, before testing different observer angles; define `d=q(t_b)-q(t_a)`.

The inherited observational slope `m=-0.308` defines `e=(1,m,0,0,0)/sqrt(1+m^2)` and `n=(-m,1,0,0,0)/sqrt(1+m^2)`. The inherited v18.52n PCA hidden direction is `h=(0,0,h3,h4,h5)` with unit norm. This `h` is a **chosen observer direction**, not a proved dynamical hinge.

For inclination alpha, set

`w_alpha=cos(alpha)*n + sin(alpha)*h`,
`d_alpha=-sin(alpha)*n + cos(alpha)*h`.

The screen is `span(e,w_alpha)`. The signed depth is `d.d_alpha` and there remain **two additional invisible 5D dimensions** orthogonal to `span(e,n,h)`.

The exact Euclidean identity is

`||d||^2 = (d.e)^2 + (d.w_alpha)^2 + (d.d_alpha)^2 + ||d_perp||^2`.

Thus `D_screen²+D_unseen²=D_5D²` for the **same** `d` regardless of angle. The velocity mismatch `||v(t_b)-v(t_a)||` is measured but does not make the relative vector an absolute orbit.

Time pair anchors: two pairs imported from v18.53d at 0° and +17.04°, frozen before the comparison; a third picked by minimizing full 5D separation in the predeclared early/late windows. For each anchored pair, five prescribed angles `(-35,-17.04,0,+17.04,+35)` are evaluated without searching again for endpoints.

If `A=d.n` and `B=d.h`, then `d.w_alpha = A*cos(alpha)+B*sin(alpha)`. An angle chosen *after seeing d* can set this component to zero, leaving `D_screen,min=|d.e|` (modulo the 180° ambiguity). This is an **observer-selection counterfactual** and not an independently identified physical angle. No false-closure claim is inferred from it unless it passes a predeclared, independent test.

## B. Native 10D to absolute 5D and the kinetic hinge candidate

Use the authenticated historical v18.49h exact kinetic embedding:

`y=(ln r1,v1,theta1,omega1,ln r2,v2,theta2,omega2,z,vz)`, `x=sigmoid(z)`, `m_x=1` (documented archived parameter).

`Q=(r1 cos(theta1),r1 sin(theta1),r2 cos(theta2),r2 sin(theta2),sqrt(m_x)*x)`;

`V=(r1(v1 cos(theta1)-omega1 sin(theta1)),r1(v1 sin(theta1)+omega1 cos(theta1)),r2(v2 cos(theta2)-omega2 sin(theta2)),r2(v2 sin(theta2)+omega2 cos(theta2)),sqrt(m_x)*x(1-x)*vz)`.

For the **declared flat kinetic Lagrangian** `T=0.5*V.V`, the kinetic Legendre moment is `P^(T)=V`. The infinitesimal `SO(2)` rotation of the **chosen** `n-h` plane has the kinetic sector momentum map

`J_hinge^(T)=(Q.n)(P^(T).h)-(Q.h)(P^(T).n)`.

In a full driven Lagrangian, nontrivial velocity-linear terms or nonconservative generalized forces may change momentum definitions and introduce torque. Neither full action invariance nor dynamic symmetry of this particular PCA observer plane has been proved. `J_hinge^(T)` is therefore a **kinetic canonical candidate**, not an automatically conserved Q550 charge.

Independently, `mathcal J= r1^2 omega1-r2^2 omega2` is the exact conserved counterrotation first integral of the **authenticated angular equations**. The audit compares absolute time series to falsify an unwarranted identity of these quantities. A finite-difference `dJ_hinge/dtau` on downsampled native states is only a numerical diagnostic, not an analytic torque certificate.

## C. Scope of folds and radius envelope

Reference files preserve the v18.53c/v18.53d authenticated sample fold radii and velocities under two κ₀ rules. However these are **equilibria at different signed drives** and do not give a time series of the same orbit. The long v18.52l native states used a historically shared κ₀. They cannot be interpreted as measuring an evolution from one dual-κ fold equilibrium to the other. The critical smaller radius is not a theorem of global dynamic noncollapse.
