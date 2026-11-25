# Mathematical protocol — Q550 v18.53j

## Historical/model scope

Native 10D states are rooted at `xi=12.8`. The authenticated source `v18.51s` supplies the Q550 RHS and historical external phase construction. Native `v18.52l` supplies fixed independent seed states; `v18.53i` supplies independent true-IVP baseline summaries through xi=51.2. The continued source clock for xi>25.6 is a **model protocol extrapolation**, not observed historical forcing. `d tau/d xi = omega/x` is checked, not inverted. No full 11D autonomous return map is claimed here.

## Residual-torque bivectors

In the physical 5D kinetic coordinates `Q` with force components from the historical RHS:

\[
 B_{\rm rad}(\tau) = Q\wedge F_{\rm radial\ potential},\qquad
 B_{\rm phase}(\tau)=Q\wedge F_{\rm phase\ coupling},\qquad
 B_{ij}=Q_iF_j-F_iQ_j.
\]

\[
 U[a,b]=\int_{\tau(a)}^{\tau(b)}B_{\rm rad}\,d\tau,\quad
 C[a,b]=\int_{\tau(a)}^{\tau(b)}B_{\rm phase}\,d\tau,\quad
 R[a,b]=U[a,b]+C[a,b].
\]

For cut points `b0<b1<b2<b3`, define `R1=R[b0,b1]`, `R2=R[b1,b2]`, `R3=R[b2,b3]` as real skew 5×5 matrices. Norms/angles use the invariant Frobenius product `<A,B>=sum_ij Aij Bij`; a common global orthogonal 5D coordinate change preserves all norms and angles. Angles near 120° are in this **10-dimensional space of skew 5×5 components**, not ordinary spatial 120° directions.

\[
 S_3=1-\frac{\|R_1+R_2+R_3\|_F}{\|R_1\|_F+\|R_2\|_F+\|R_3\|_F}.
\]

`S3≈1` describes a small net relative to separately accumulated signed bivectors. It is neither energy conservation, nor closed 5D trajectory, nor the conserved native charge `J=r1² ω1-r2² ω2`. The independently reported radial-phase score `Srad-phase=1-||sum(Uj+Cj)||/(||sum Uj||+||sum Cj||)` differs from `S3` because its denominator is different.

## Critical additive null hypothesis

At **fixed external start/end**, the sum of the three integrals is exactly

\[
 R_1+R_2+R_3=R[b_0,b_3],
\]

whatever internal cuts `b1,b2` are chosen. Therefore the absolute net norm is a mathematically guaranteed invariant of internal-boundary jitter. Our script actively tests this equality as a **software identity**, not a physical result. The denominator of `S3`, individual norms, pairwise 120° angles, singular strengths and principal-plane polygon diagnostics may still change with cuts. Varying the **total origin** and/or **total span** is necessary to attempt falsification of the net near-closure itself.

## Principal planes / discrete geometric return

For each real skew tensor `Rj` in 5D, a real Schur decomposition yields up to two singular 2-planes with strengths `sigma1 >= sigma2 >= 0` and a 1D kernel. Individual axes are gauge-dependent. A principal-plane projector `P=UU^T` is gauge-invariant except near singular-strength collisions. We report a relative singular gap and undefined return angle when a principal plane is not identifiable.

For one selected mode, take frames `(U1,U2,U3)` of the rank-2 planes, and close the polygon **by convention** with `U3 -> U1` along the shortest principal-angle geodesic in the Grassmannian. The basis of the next plane is aligned to the transported previous frame using the orthogonal Procrustes polar factor of `U_next.T @ U_prev`. The resulting two-by-two return `W=U1.T @ U1_return` is orthogonal. `acos(trace(W)/2)` gives an unsigned angle for `det(W)>0`; if the minimum singular overlap of adjacent frames is too small, the return is declared undefined. This O(2)-frame-gauge-invariant **constructed geodesic polygon** measures the geometry of the three observed *principal-plane subspaces*. It does NOT supply a physically generated connection or a genuine temporal holonomy of a closed orbit; the physical evolution is open. Do not call this a material-memory certificate.

## Pre-registration and numerical checks

- Three actual κ rules: shared, original and autonomous; two original native historical seed families.
- Whole triplets: origins `[12.8,13.6,14.4,15.2]`, each of equal widths `[9.6,10.8,12.0,12.8,13.6]`.
- Fixed-endpoint jitter: `b0=12.8`, `b3=51.2`, `b1=25.6+delta1`, `b2=38.4+delta2`, `delta1,delta2` in `[-2,-1,0,1,2]`.
- Numerically accepted near-closure `S3>=.90`, near equiangular max angle error <=8 degrees, and nearly equal magnitudes CV<=.10. Thresholds are classifications, not proof of a law.
- Cumulative trapezoid in **material time**, linear interpolation of cumulative integrals to predetermined off-grid edges, explicit additivity check, and independent direct interpolated-edge trapezoid audit for off-grid subintervals.
- Source RHS error <2e-10, native replay <2e-4, charge drift <5e-9, refined shared-κ IVPs agree <3e-7, exact genealogy and SHA checks of all embedded scientific ancestors.

## Interpretation hierarchy

1. *Numerical identity:* fixed-endpoint additivity and its rotation invariance.
2. *Authenticated observed diagnostic:* 120-degree-like bivector-space triplet and near-closure under a specified clock/interval/κ, provided all replay/source checks pass.
3. *Window robustness:* empirical persistence across predefined origins, durations, and internal cuts; may fail without invalidating underlying historical physics.
4. *Constructed observer-space polygon transport:* singular-plane geodesic polygon return diagnostic, NOT a physical loop holonomy.
5. *Unproven:* full-ramp triple rhythm, material return to identical state, new Noether charge, global attractor or physical three-step cycle.
