# Mathematical protocol — v18.52x

## Physical and external-clock state

The physical state is `y=(l1,v1,theta1,omega1,l2,v2,theta2,omega2,z,vz)` with `r_i=exp(li)`. Native v18.52l column 11 is accumulated material time `tau`, **not external phase** `Theta`. The authentic external phase is recovered from the historical v18.51s source `ctx['theta']`.

The conserved angular difference is

\[
J(y)=e^{2l_1}\omega_1-e^{2l_2}\omega_2, \qquad \dot J=0.
\]

The continuous symmetry is simultaneous counterrotation `theta1 -> theta1+epsilon`, `theta2 -> theta2-epsilon`. The quotient coordinate `alpha=theta1+theta2 (mod 2pi)` is gauge invariant. Fixing `J=J_*` yields an eight-dimensional chart

\[
u=(l_1,v_1,\omega_1,l_2,v_2,z,v_z,\alpha), \qquad
\omega_2={e^{2l_1}\omega_1-J_*\over e^{2l_2}}.
\]

A pure gauge angle difference is removed **only** from the configuration; both angular velocities remain physically relevant before fixing the charge. Fixing charge eliminates `omega2` from the chart by the exact constraint, not by an arbitrary projection. Euclidean distances on these mixed-unit chart coordinates are descriptive, **not** a kinetic metric.

## External-phase Poincaré section

On the long return branch the authenticated signed clock is `s=-1` and the original source gives a strictly decreasing `Theta(xi)`. The physical source generator is driven by `Theta mod 2pi`; therefore define a section at a fixed value of `Theta mod 2pi`, with the **negative crossing direction**. Match endpoints with `Theta(xi_b)=Theta(xi_a)-2pi`, not merely a configuration phase match. The pairs

\[(u_a,J_a)\mapsto (u_b,J_b),\qquad J_b\simeq J_a\]

are authentic *one-step observations*. The finite native interval has fewer than two complete phase windings; its observations do not provide iterates `P(P(u))` at the same section. Using 17 distinct anchor phases yields 17 **different section phase choices**, not 17 sequential iterations of a single map.

The reduced differences use `wrap(alpha_b-alpha_a)`. Each of the other seven components is compared directly. A small sum-angle phase mismatch alone does not certify a periodic orbit if hidden velocities change.

## Charge stratification / fixed-J identification

Different historical kicks are **not on the same conserved-J leaf as the reference**. Each trajectory may be used for its own valid one-step return pair, but cross-trajectory differences cannot identify `D P_{J_ref}`. The four signed native kicks span at most four directions in an eight-dimensional charge-fixed chart and in practice predominantly one nominal kick family. Base/refined copies are the same perturbations for numerical convergence, **not extra independent directions**. A singular-value rank is reported diagnostically but no full transverse eigenvalues are inferred.

The central kick secant

\[d_\epsilon u_{in/out}={u(+\epsilon)-u(-\epsilon)\over2\epsilon}\]

and its norm ratio are labeled **mixed-J family responses**. They need not equal a fixed-J tangent derivative, because `dJ/d epsilon != 0`.

## Hierarchy of conclusions

1. One-step phase-gated full-state return measured: empirical section map observations only.
2. Exact relative periodic orbit: requires vanishing full reduced mismatch with controlled numerical error on a genuine repeat; near-zero is not exact.
3. Fixed-J 8D derivative: requires >=8 independent *same-J* tangent directions (or a rigorously implemented analytic variational equation) and rank/conditioning checks.
4. Attracting relative limit cycle: additionally requires appropriate transverse multipliers inside the unit circle, with neutral symmetry/charge directions properly removed, and an actual periodic orbit.
5. Full-system conservative dynamics/global energy conservation: **not** implied by the Noether charge, a phase match, or a limit cycle on a dissipative forced leaf.
