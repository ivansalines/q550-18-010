# Q550 v18.52i — Mathematical protocol

The matched-material-time relative trajectory in the physically justified 5D
Cartesian embedding from v18.49h is

\[q(\tau)=Q_{\mathrm{kick}}(\tau)-Q_{\mathrm{reference}}(\tau),\qquad
v(\tau)=\dot q(\tau),\qquad a(\tau)=\ddot q(\tau).\]

The physical kinematic identities remain

\[D=\tfrac12q\cdot q,\qquad G=D'=q\cdot v,\qquad
C=D''=\|v\|^2+q\cdot a=P-I,\quad P=\|v\|^2,\ I=-q\cdot a.\]

**A new local minimum** is inferred only from a sampled sign-changing bracket
with `G: negative -> positive`, and requires refining that zero on the authentic
nonlinear dense flows. All tangencies not bracketed by the finite grid remain
uncertified. New `C` sign changes alone are *not* minima.

Geometric decomposition (only where `r=||q||` is above a reported floor):

\[e=q/r,\quad \dot r=e\cdot v=G/r,\quad
v_{\perp}=v-\dot r e,\quad \Omega=\|v_\perp\|/r.\]

Compute the angle
\[\alpha(\tau)=\arccos(e(\tau)\cdot e(\tau_h)),\]
where `tau_h` is the **first common material time available after the h seam**,
not an assertion of equal global xi across kicked and reference clocks.
Also store cumulative directional path \(\int\Omega\,d\tau\) and the actual
5D difference vectors, so angles do not depend on a plotted 2D projection.

One can have `D'=0` but `v_perp != 0`: an apparent radial turning point does
not make the 5D vector stationary. Likewise `D(t1)=D(t2)` does not imply
`q(t1)=q(t2)`. These are the primary closure-discrimination tests.

The generator must use the exact global xi argument for the drive and preserve
the extra material clock state. Refined independent 10D IVPs check observed
bracketed roots and physical vectors. The output does not claim an exhaustive
continuous-time proof or physically exact periodicity.
