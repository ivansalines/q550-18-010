# Mathematical protocol — Q550 v18.52y

## 1. Exact Noether-charge leaf and quotient

The authenticated 10D physical state is `Y=(l1,v1,t1,w1,l2,v2,t2,w2,z,vz)` with positive radii `ri=exp(li)`. The exact first integral is `J=r1^2*w1-r2^2*w2`. The counterrotation symmetry is `t1->t1+eps`, `t2->t2-eps`, with invariant `alpha=t1+t2`. Fix `J=J0` and use the eight-chart `u=(l1,v1,w1,l2,v2,z,vz,alpha mod 2pi)`. Choose any fixed gauge `g=t1-t2`. Then

`L(u;J0,g) = (l1,v1,(alpha+g)/2,w1,l2,v2,(alpha-g)/2,(exp(2*l1)*w1-J0)/exp(2*l2),z,vz)`.

Its columns are linearly independent and tangent to `J=J0`. In particular

`dw2/dl1=2 exp(2*(l1-l2))*w1`, `dw2/dw1=exp(2*(l1-l2))`, `dw2/dl2=-2*w2`, `dt1/dalpha=dt2/dalpha=1/2`.

All eight columns are verified with finite differences and `∇J·DL=0`. This imposes J *before propagating*, not after.

## 2. True stroboscopic map of the late return protocol

On the authenticated late return branch `s(t)=cubic_clock_sign(t)=-1` and `Theta(xi)=ctx.theta(ctx.event+xi)`, so `Theta(xi+T)-Theta(xi)=-2*pi` with `T=2*pi/engine.omega`. Integrate

`dY/dxi=F_phys([Y,Theta(xi)], kappa_fixed, source_params, source_slope, s=-1, source_support, source_audit)`,

using the authenticated `rhs_from_original`, not a guessed or linearized model, and the 10D native initial state at `xi=12.8`. Define `P_J(u)=reduce8(Flow_T(L(u;J,g)))`. The observed native endpoint at `xi=12.8+T` is an independent replay check. Continue to `xi=12.8+2*T` for the second genuine stroboscopic image; compare `P_J(P_J(u))` with the direct two-period flow to verify late periodicity.

The map is a source- and protocol-conditional mathematical object, not a claim that the historical complete forward–return drive was globally periodic.

## 3. Transverse derivative

For each of eight coordinate directions `e_i`, set `ui±=u±h_i e_i` and lift both to the same `J0`. The numerical Jacobian is the centered difference

`DP_J[:,i] = wrap_difference(P_J(ui+),P_J(ui-)) / (2*h_i)`.

Use `h_i=eps*max(1,abs(ui))` for `eps=1e-4` and `4e-5`; compare the matrices. Plot coordinate singular values as **coordinate-dependent diagnostics only**. Eigenvalues of a Jacobian at an *unclosed* state are NOT Floquet multipliers and MUST NOT be interpreted as stability of a periodic orbit.

## 4. Relative closure and acceptance

Seek `R(u)=wrap_difference(u,P_J(u))=0` using a bounded damped Newton update with cap on weighted step and backtracking. A root is accepted only if independent recomputation finds `||R(u*)||_8 < 1e-7` and `||wrap_difference(u*,P_J^2(u*))||_8 / max(1,||u*||) < 2e-7`, together with checked composition and charge. A closure is modulo counterrotation and the 2π angular periodicity. A successful root still requires two-epsilon Jacobian agreement before reporting conditional Floquet multipliers; transverse attraction is only stated when all moduli are strictly below one by a numerical margin.

## 5. Forbidden inferences

- Do not reuse mixed-charge kicks for the fixed-charge Jacobian.
- Do not label a non-fixed-point Jacobian spectrum as Floquet.
- Do not infer global dissipation balance or total-energy conservation from J.
- Do not infer the complete forward–return cycle from the late constant negative-drive stroboscopic map.
- Do not treat the chart norm as a physically invariant metric.
