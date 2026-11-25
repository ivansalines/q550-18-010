# Q550 v18.53f — Mathematical protocol

All derivatives in this note are with respect to **physical material time** tau, not the auxiliary xi. Historical v18.52l stores tau as its **eleventh state column**. The phase Theta is separately provided from authenticated v18.51s context `theta(event+xi)`. The native 10D dynamics use the **shared** `kappa_fixed=0.1666331385421536`; historical original/autonomous kappa rules are not independent native trajectory parameters here.

## 1. Flat kinetic embedding

For material state `y=(l1,v1,theta1,w1,l2,v2,theta2,w2,z,vz)`, let `r_i=exp(l_i)`, `x=sigmoid(z)`, `xd=x(1-x)vz`, `mx=1` in the authenticated source.

`Q=(r1 cos theta1,r1 sin theta1,r2 cos theta2,r2 sin theta2,sqrt(mx)*x)` and `V=(r1(v1 er1+w1 et1), r2(v2 er2+w2 et2), sqrt(mx)*xd)`. For the declared flat kinetic sector, `P^(T)=V`.

For an ordered fixed orthonormal plane `(n,h)`, define

`J_hinge^(T)=(Q.n)(V.h)-(Q.h)(V.n)`.

This is a *kinetic-sector* moment. It is not proven to be the full Q550 canonical Noether generator. The separate exact charge is `J_rot=r1^2*w1-r2^2*w2`.

## 2. Exact fixed-plane torque identity

At fixed `(n,h)`, with `A=dV/dtau`,

`dJ_hinge^(T)/dtau=(Q.n)(A.h)-(Q.h)(A.n)`.

The `V wedge V` term vanishes identically. The physical 5D acceleration is decomposed into literal model-source force channels:

1. Radial potential: `F_U,i=-U'(ri) er_i`, with `U'(r)=2r-2.4r^2+0.8r^3-e_reserve^2/r^3`.
2. Phase coupling (radial and tangential): `F_c,i=-kappa*x*rj*cos(Delta) er_i + kappa*x*rj*sin(Delta) et_i`, `Delta=theta1+theta2-phi0-Theta_external`.
3. Radial damping: `F_d,i=-gamma_rho*ri*vi er_i`.
4. Logistic auxiliary barrier in 5D: `F_B,5=-B*(2x-1)/(sqrt(mx)*x*(1-x))`.
5. Auxiliary active source in 5D: `F_mu,5=mu_x*r1*r2*sin(Delta)/sqrt(mx)`.
6. Auxiliary damping in 5D: `F_dx,5=-gamma_x*xd/sqrt(mx)`.

`A=sum(F_channel)` must agree with acceleration computed by differentiating the 5D embedding along the staged **original 11D source RHS** after converting `dy/dxi` to `dy/dtau` by `x/omega`. The last source equation is checked against `Theta_xi=omega*cubic_clock_sign`. The same signed clock controls the historical external phase. A direct signed-drive *force* is **not** postulated: at a frozen state and frozen external phase the first ten RHS components do not change upon reversing the clock-sign argument; the clock affects *subsequent* force evolution through Delta.

Linearity of the wedge gives six separate source torque channels and their sum. Source decomposition is conditional on the documented chart and parameterization; it is not a unique decomposition of physical interaction energy.

## 3. Moving observer and internal basis rotations

For time-varying orthonormal `(n(tau),h(tau))`:

`dJ_obs/dtau = wedge(Q,A;n,h) + wedge(Q,V;ndot,h) + wedge(Q,V;n,hdot)`.

The final two terms are **observer transport**, not physical torque. For a basis rotation *inside the same oriented plane*, `n_a=cos(a)n+sin(a)h`, `h_a=-sin(a)n+cos(a)h`, the two geometric terms cancel and `J_obs=J_fixed` exactly. As a controlled test of an observer plane that genuinely changes in R5, choose an additional `e` perpendicular to `(n,h)`, and tilt `n_beta=cos(beta)n+sin(beta)e`, `h_beta=h` with a fixed, predeclared sinusoidal schedule. This is a probe of observer dependence, NOT a claim that the physical hinge itself is rotating.

## 4. Hidden advance at identical endpoint times

For the archived relative vector `dQ(tau)=Q_kick-Q_ref` and `dV=d(dQ)/dtau`, for a *fixed* observer screen `(e,w_a)` and hidden direction `d_a`, the endpoint full squared distance splits exactly into screen squared plus hidden squared. Verify independently that `dQ(t_b)-dQ(t_a)=integral_(a,b) dV dtau`, and project that identity onto each predeclared hidden direction. Never relabel this as an absolute Q550 orbital closure.

For a moving hidden basis vector `d(tau)`, the observed signed hidden-coordinate change satisfies

`(dQ_b.d_b)-(dQ_a.d_a)=integral[dV.d]dtau + integral[dQ.ddot]dtau`.

First term: actual relative motion resolved along the moving direction. Second: motion *of the observer frame*. The same physical endpoints remain fixed for every predeclared view; no angle is optimized after seeing the data in this audit.

## 5. Decision boundary

- Full pass: original source and clock, fixed native kappa, force identity, fixed/moving torque ledgers, in-plane invariance, authenticated charge, and relative hidden-advance integral all pass.
- Partial: historical parent missing; report native kinetic and projection identities only, no authenticated force-source attribution.
- Reject: any SHA mismatch, source parameter drift, force identity fail, or quantitative ledger fail.
- NEVER conclude full hinge Noether conservation, new global orbit closure, or physical transport of a unique hinge basis from these data alone.
