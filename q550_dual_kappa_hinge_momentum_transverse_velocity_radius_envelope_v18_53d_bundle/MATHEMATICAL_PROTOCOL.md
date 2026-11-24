# Q550 v18.53d — mathematical protocol

## A. Authentic dual-kappa rotor mechanics

Inputs: the SHA-pinned authenticated output of v18.53c, not recreated dynamics. Along each `Delta=pi` equilibrium fold at fixed charge `J=-1.441016445824766`, use the recorded `r1,r2,omega1,omega2` to evaluate:

\[
 v_{\theta,i}=r_i\omega_i,\qquad J_i=r_i^2\omega_i,\qquad \mathcal J=J_1-J_2.
\]

The `min(r1,r2)` at a fold is the smaller *equilibrium* radius. Different kappa values have different signed drive values at their folds. The original/autonomous kappa comparison is not an evolution at fixed drive; a separate same-drive `s=-1` radius comparison is inherited from v18.53c. No continuous interval bound or dynamical noncollapse theorem is asserted by this audit.

## B. One authentic relative 5D trajectory, multiple observers

Input: the SHA-pinned v18.52l `extended_horizon_5d_vectors_and_traces.npz` containing matched-time `(dQ,dV)` for eight perturbation-minus-reference traces. The observed slope `m=-0.308` and transverse hidden PCA vector come from the v18.52n **fixed-hypothesis observer**; the visually observed four-point slope has **not** been independently fitted.

Define the orthonormal hinge axis, planar normal, hidden reference direction:

\[
 e=\frac{(1,m,0,0,0)}{\sqrt{1+m^2}},\quad
 n=\frac{(-m,1,0,0,0)}{\sqrt{1+m^2}},\quad
 h=(0,0,h_3,h_4,h_5).
\]

For observer inclination `alpha`, define the image plane span`{e,w_alpha}` with

\[
 w_\alpha=\cos\alpha\,n+\sin\alpha\,h,\qquad
 d_\alpha=-\sin\alpha\,n+\cos\alpha\,h.
\]

For two different material times, write `Delta q = dQ(t_b)-dQ(t_a)`. Compute:

\[
 D_{2}^2=(\Delta q\cdot e)^2+(\Delta q\cdot w_\alpha)^2,\quad
 D_{5}=\|\Delta q\|,\quad
 D_{\rm unseen}^2=D_5^2-D_2^2.
\]

The signed *one-dimensional* viewing-depth displacement is `Delta q . d_alpha`. Other hidden orthogonal components may contribute to `D_unseen`. In particular, `D2=0` never implies `D5=0`. The code selects the smallest 2D distance between predeclared early/late windows, at each **predeclared** observer angle (-35,-17.04,0,+17.04,+35) — it does not tune the angle to a selected pair. Thresholds define only a descriptive candidate for apparent closure. A velocity mismatch is checked separately. Two samples on different times do not establish orbit recurrence or periodicity.

## C. Hinge kinetic wedge versus canonical moment

For the **relative** vectors `(dQ,dV)`, compute a signed *kinematic area rate*

\[
 L_{nh}^{\rm kin}=(dQ\!\cdot n)(dV\!\cdot h)-(dQ\!\cdot h)(dV\!\cdot n),\quad
 R_{nh}=\sqrt{(dQ\!\cdot n)^2+(dQ\!\cdot h)^2},\quad
 v_{\theta,nh}=L_{nh}^{\rm kin}/R_{nh}.
\]

The wedge is invariant under passive SO(2) frame rotations *within its own n-h plane*, but it is not an observer-invariant 2-form over arbitrary plane choices. The 2D screen wedge `L_{e,w_alpha}` is inclination dependent. The physically canonical candidate,

\[
 J_{\rm hinge}=(q\cdot n)(p\cdot h)-(q\cdot h)(p\cdot n),
\]

**cannot be evaluated** without the genuine absolute 5D `(q,p)`, the model's kinetic/canonical momentum map and a validated hinge plane. **Do not replace `p` with `dV` silently**. In particular, neither 5D kinematic wedge is the proven conserved charge `J1-J2` by assumption.

## Decision controls

- Verify copied v18.53c report / fold CSV / fold states against its authenticated output manifest.
- Verify v18.52n observer report against authenticated output manifest.
- Verify v18.52l 5D data against SHA-256 `d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`.
- Require each fold point to have authentic RHS `native_checked` and `native_source_match`.
- Distinguish the *two model kappa values* of the equilibrium fold from the *single shared kappa* used to generate the historical v18.52l 5D traces.
- Never certify closure, canonical hinge moment, topological winding or a global minimum radius from an observer projection.
