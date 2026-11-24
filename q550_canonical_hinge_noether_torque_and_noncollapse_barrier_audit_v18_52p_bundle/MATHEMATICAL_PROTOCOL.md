# Mathematical protocol — Q550 v18.52p

## 1. Source and coordinate geometry

The authenticated v18.49h kinetic embedding defines (for each genuine 10D state)

`Y=(ell1,v1,theta1,w1,ell2,v2,theta2,w2,z,vz)`, `r_i=exp(ell_i)`, `x=sigmoid(z)` and

`Q=(r1 cos(theta1), r1 sin(theta1), r2 cos(theta2), r2 sin(theta2), sqrt(m_x)*x)`,

`V=(r1(v1 e_r1+w1 e_t1), r2(v2 e_r2+w2 e_t2), sqrt(m_x)*x(1-x)*vz)`.

The natural kinetic energy is `T=0.5*||V||_5²`; in a **natural Cartesian Lagrangian with no velocity-dependent interactions**, `P=V`. This is an explicit **kinetic Legendre proxy**, not proof that the driven Q550 flow has such a canonical action. For the first and second physical polar pairs, the kinetic proxy is `J1=Q1 V2-Q2 V1=r1²*w1` and `J2=Q3 V4-Q4 V3=r2²*w2`. The first two components are the physical Cartesian coordinates; the 5D relative separation `q=Qkick-Qreference` is **not an independent physical body**.

## 2. Canonical rotation and nodal-axis sign

With observed slope `m=-0.308`, choose `theta=atan(m)<0` and passive rotation

`R(theta)=[[cos(theta),sin(theta)],[-sin(theta),cos(theta)]]`.

Then the first row equals the nodal unit vector `(1,m)/sqrt(1+m²)`, and the second row is the planar normal `(-m,1)/sqrt(1+m²)`. Apply *the same R* to both `Q12` and kinetic `P12`. The 4×4 `S=diag(R,R)` satisfies `S^T Ω S=Ω` and exactly

`Q1' P2' - Q2' P1' = Q1 P2 - Q2 P1`.

This equality under coordinate change is **not** `dJ/dtau=0`. Also, the user's +17.04° passive display matrix aligns an axis of *positive* slope. For `m<0`, the basis rotation must use `theta=atan(m)` with negative sign (or explicitly reconcile active versus passive conventions). Paper tilt ±17.04° into hidden Q3–Q5 observer directions is a distinct operation from the canonical planar rotation.

## 3. Actual Noether requirement

With canonical `(Q,P)`, a planar active infinitesimal rotation has `δQ=(-Q2,Q1)ε`, `δP=(-P2,P1)ε`. Its generator is `J=Q1 P2−Q2 P1`. For a truly canonical Hamiltonian flow,

`dJ/dtau = {J,H} = Q2*(∂H/∂Q1) - Q1*(∂H/∂Q2) + P2*(∂H/∂P1) - P1*(∂H/∂P2)`.

For a natural Cartesian isotropic kinetic `H=T+V(Q,...)` the kinetic terms cancel, leaving the planar potential torque `Q2 V_Q1 - Q1 V_Q2`. Add generalized nonconservative torque if the original model includes forcing, dissipation, noncanonical clocks or velocity-dependent terms. An unmodified, actual action/RHS equivariance check is required for a Noether conservation certificate. **This audit deliberately does not invent `H`.**

## 4. Independent on-trajectory kinetic balance

The relative physical 5D source provides `q=dQ`, `v=dV`, `a=da` on the same material clock. For planar components `p=(q1,q2)`, set

`r=||p||`, `Jkin=p1*v2-p2*v1`, `torque=p1*a2-p2*a1`.

Then the exact product rule states

`(d/dtau) Jkin = torque`.

The script compares a numerical derivative of Jkin to the archived acceleration torque, to independently check data consistency. It also reconstructs *absolute* reference/kick kinetic angular quantities from the authentic 10D polar states and aligns `Jkick` and `Jreference` at **common material time τ** via monotone interpolation. For two distinct flows,

`Jkick-Jreference != (Qkick-Qreference) x (Pkick-Preference)` in general;

the mismatch is a measured feature, not a numerical error. A full two-flow canonical relative transformation, even if it exists, has `p_rel=(Pkick-Preference)/2` for a pair of equal-mass independent bodies, not `p_rel=Pkick-Preference`, and it does not automatically generate a conserved relative J.

## 5. Planar decomposition, conditional lower bound and 5D bivector

For `r>0`, define `rdot=(p·v2)/r`. The exact identities are

`0.5*||v2||² = 0.5*rdot² + Jkin²/(2r²)`

`rddot = Jkin²/r³ + (p·a2)/r`.

The first identity motivates a **candidate** centrifugal barrier. It is NOT an effective potential with conserved `J` unless the force law and kinetic structure justify it. Cauchy-Schwarz implies `r>=|Jkin|/||v2||` at each time with nonzero velocity. Using `jmin=min_grid|Jkin|` and `vmax=max_grid||v2||` yields an **a posteriori bound only on the samples**. A genuine no-collapse theorem would need a continuous/global lower bound on J and an upper bound on kinetic speed or energy; sampled extrema alone cannot provide either for future times.

The full 5D angular bivector has 10 components `Lij=qi*vj-qj*vi` and `dLij/dtau=qi*aj-qj*ai`. Its norm satisfies `||L||²=||q||²||v||²−(q·v)²`. A conserved 5D bivector would require an appropriately rotationally symmetric 5D system, which is NOT assumed here.

## 6. Visual phase gates and open-path winding

For the historical 5D observer axis `h` pinned by the v18.52o report, let `n=(-m,1,0,0,0)/sqrt(1+m²)` and

`w±(tau)=cos(17.04°)*(n·q) ± sin(17.04°)*(h·q_hidden)`.

Zero crossings are extracted without fitting to angular momentum; signed `J`, `r2`, `r5` and radial balance are evaluated at the same gates. `turns=[unwrap(atan2(q2,q1))]_end-start/(2π)` is valid for an **open oriented segment** if `r2>0`; it is not an integer winding/topological charge without an admissible closed loop and a nonvanishing radius during a homotopy.

## 7. Numerical controls

- SHA-256 ancestry for 5D and native 10D sources; exact inherited v18.49h kinetic file.
- Tests: exact canonical J under rotation, symplectic matrix, natural central conservative toy (J constant and radial cancellation), noncentral torque toy, 5D bivector identity, native polar kinetic identity and zero-J collapse-allowing toy.
- Do not confuse `tau` with event-clock `xi` or align independent solutions by sample indices when comparing kick/reference quantities; use the native physical material clock.
- `sample_stride>1` always labels a report `DEVELOPMENT_DIAGNOSTIC_ONLY`; normal server execution uses full resolution and requires native 10D data.

## 8. Chart versus physics (no false noncollapse certificate)

The exact original chart expresses `r_i=exp(ell_i)>0` for finite `ell_i`. This coordinate identity is not a theorem that `r_i` stays uniformly away from zero; approaching `ell_i -> -infinity` could still correspond to approach to a physical singular radius. Similarly `q=Qkick-Qreference` can approach zero even when each absolute radius remains positive. Physical noncollapse must always specify whether it refers to a native rotor radius, a nodal plane radius of *relative separation*, or the full 5D separation norm.
