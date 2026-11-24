# Decision rules — v18.52e

1. **Numerical prerequisites:** authenticated v18.52d sendback and exact v18.51s ancestry; 8 kick IVPs; no solver failures; no nonfinite values; physical clock strictly increasing; baseline reproduction better than `1e-6` in kinetic coordinates; dot-product identities and clock synchronization within documented tolerances.
2. **Relative speed:** compute the signed power `P_rel=dV·da`, its negative fraction, integrated power, and ratio of terminal to initial `|dV|`. Negative power does not imply physical dissipation.
3. **Position:** compute signed separation rate `G_rel=dQ·dV`, its sign and the terminal distance. Distinguish position restoration from velocity-difference contraction.
4. **Lateral turning:** report `|a_perp|/|V|²`, change of curvature and the acceleration difference's lateral fraction. Curvature change may be of either sign.
5. **Initial-point handling:** at t=0, dQ=0. Do not classify `dQ·da` as an inward restoring effect there. Report the qualified sample count and thresholds.
6. **No global claims:** a finite-window linearized or nonlinear kick diagnostic at one parameter and chart does not prove material memory, asymptotic stability, physical dissipation, or charge-conserving realizability.
7. **Reported numerical uncertainty:** spline differentiation checks are numerical comparisons, not rigorous interval envelopes; all sample-based sign classifications remain qualified.
