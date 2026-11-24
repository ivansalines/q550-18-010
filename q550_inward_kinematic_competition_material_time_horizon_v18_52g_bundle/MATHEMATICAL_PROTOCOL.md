# v18.52g mathematical protocol

Let both trajectories have the same historical material parameter κmid, the same original nonlinear 10D field and the same seed; only the initial finite lateral Cartesian velocity kick differs. `xi` is the external-event-relative integration variable, and the material clock τ is integrated as an eleventh state with `dτ/dxi=c(Y,xi)>0` based on the historical physical clock formula.

Each branch is evaluated at its own `xi_branch(τ)` by root-finding on the numerical dense 11D solution. This ensures identical material time, though the external drive phase xi may differ; the audit records that mismatch.

In the genuine physical Cartesian 5-vector embedding, define

```
q = Qkick(τ)-Qreference(τ)
v = Vkick(τ)-Vreference(τ)
a = akick(τ)-areference(τ)
D = 0.5 q·q
G = q·v = dD/dτ
P = v·v >= 0
I = -q·a
C = P-I = dG/dτ = d²D/dτ²
R = I/P  (only for P>0)
```

Threshold hierarchy:

- `I>0`: projection of relative acceleration points inward; *not* sufficient for slowdown or return.
- `C=0` (or `R=1`, when P>0): convexity changes only if a sign crossing occurs; may still have `G>0`.
- `G=0`: stationary squared separation; a positive→negative crossing gives a local maximum of D and a subsequent decrease.
- Subsequent smaller `||q||`: finite-window position approach, not necessary 10D state recovery or long-time stability.

Use DOP853 dense solutions and Brent roots on matched τ; do **not** locate events exclusively from a PCHIP/spline of the old f/e arrays. Verify `xi<=0.2` against parent history, local identities and a refined independent integration. Root search is limited by the material-time grid, which is **not an interval enclosure**. `--xi-end 0.8` is an exploratory fourfold extension of the previously authenticated event-clock horizon, not a prior theoretical prediction.
