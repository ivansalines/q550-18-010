# Mathematical protocol (source-derived model equations and test hypotheses)

Source state: `y=(ln r1, v1, theta1, omega1, ln r2, v2, theta2, omega2, z, vz, Theta_ext)`, `x=sigmoid(z)`, `Delta=theta1+theta2-phi0-Theta_ext`.

The original *physical-time* equations reconstructed from `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py` include:

- `dr_i/dtau = r_i*v_i`.
- `d^2 r1/dtau^2 = r1*omega1^2 - U'(r1) - k*x*r2*cos Delta - gamma_rho*r1*v1` (and exchange indices for rotor 2).
- `d omega1/dtau = k*x*(r2/r1)*sin Delta - 2*v1*omega1`, same with 1 <-> 2.
- `d^2 x/dtau^2 = [-Vbar'(x)+mu_x*r1*r2*sin Delta -gamma_x*dx/dtau]/m_x`.
- `J=r1^2 omega1-r2^2 omega2` remains invariant identically, but total E need not be.
- `U'(r)=2r-2.4r^2+0.8r^3-e_reserve^2/r^3`, `Vbar'(x)=B*(2x-1)/(x*(1-x))`.

The historical native source coordinate xi obeys `d/dxi=(omega/x)*d/dtau` for the physical 10D state. In the late fixed-negative return protocol: `dTheta_ext/dxi=-omega`, hence `dTheta_ext/dtau=-x`. Thus a locked relative equilibrium with constant physical coordinates except the angular phases must satisfy:

`v1=v2=vz=0`, `x=1/2` on the chosen sin(Delta)=0 branch, `Delta=pi (mod 2*pi)`, `omega1+omega2=-1/2`,

`0=omega1^2 - U'(r1)/r1 + k*(1/2)*r2/r1`,
`0=omega2^2 - U'(r2)/r2 + k*(1/2)*r1/r2`,
`omega2=(r1^2*omega1-J)/r2^2`.

This is an independently rooted algebraic system of three equations for `(logr1,logr2,omega1)` at fixed `J`. Verify native 11D RHS, including the auxiliary dynamics, separately. It is not guaranteed that an 8D map fixed point is such an equilibrium; that is exactly the falsifiable question.

Use the intrinsic **autonomous reduced coordinates** `u=(logr1,v1,omega1,logr2,v2,z,vz,Delta)` after removing the external phase. For late constant-sign drive, `du/dxi = F_red(u)`, with the last component `F[theta1]+F[theta2]-F[Theta_ext]` and `omega2` determined by J. At a stationary root, the local monodromy eigenvalues are `exp(eig(D F_red)*T_xi)`. Compare these independently with the v18.52z return-map spectrum using an optimal complex pairing. This test is conditional on the parameter, leaf, phase protocol and quotient chart being *identical*.

A local phase-locked relative equilibrium is not an analytic global proof of nonlinear attraction, not a conservative Hamiltonian cycle, and does not prove the complete historical ramp forward-return is periodic.
