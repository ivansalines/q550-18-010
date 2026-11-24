# Q550 v18.52u — Mathematical protocol

## Two different operations, not one
Original autonomous source field uses internal `Theta_engine` and `dTheta_engine/dtau=+x` with coupling phase `Delta=theta1+theta2-phi0-Theta_engine`.

**Passive relabeling**: define `Theta_ext=-Theta_engine` while all first ten state coordinates are unchanged. The push-forward field necessarily has `dTheta_ext/dtau=-x` **AND** `Delta=theta1+theta2-phi0+Theta_ext`. This is an exact coordinate identity. It is not experimental evidence that the archived external clock implements this relabeling.

**Externally driven archives**: previous v18.52t uses `Theta_ext` supplied by the authenticated `ctx['theta']`, shared between kicks, and uses `Delta=theta1+theta2-phi0-Theta_ext`. The saved 10D physical dynamics are not automatically equivalent to the pushed-forward original 11D autonomous field. Test the force equations, not just `dTheta/dtau`.

The natural candidate conjugacy identity for a fixed matrix C is `F_archive(C z) = D C(z) F_engine(z)`. Audit identity and angle-sign-flip candidates on off-orbit states; do not certify full conjugacy merely from a fitting scalar or from phase-derivative agreement. A proof for arbitrary nonlinear transformations requires further derivation.

## Genealogical test
SHA-check the actual v18.51s parent, stage its ancestral engine with the exact historical loader, inspect the source implementation of `make_context` and its theta callable and the original `structural_jet`, then evaluate the external theta clock at authenticated event+xi. Source snippets and hashes are exported. No invented interpolation or reconstruction of unprovided raw sources.

## Native force discrimination
On reference series compare `d omega_i/d tau` (computed from saved states) against both signs of the coupling `Delta` using the same authenticated external phase: (a) archived `-Theta_ext`; (b) passive re-label `+Theta_ext`. The native finite-difference derivative has sampling error, report relative RMS without assuming exact zero. Neither comparison proves that the full 11D dynamics are conjugate, but it can falsify a naive sign-only coordinate interpretation.

The external reference clock shared with kicks is not an independent state governed by `Theta_dot=-x_kick`; `x_kick` differs from `x_reference`. This is a non-autonomous inherited forcing, not an 11D autonomous model with per-kick `Theta_dot=-x`.

The Noether-d'Alembert conserved charge `J1-J2` is unaffected by the phase sign since both angular momenta receive equal torque; it does not imply a uniform-in-time radial or 5D-separation noncollapse bound.
