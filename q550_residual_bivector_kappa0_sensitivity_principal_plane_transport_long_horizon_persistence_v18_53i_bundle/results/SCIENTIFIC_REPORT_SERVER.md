# Q550 — v18.53i: κ₀ sensitivity, principal-plane transport & long-horizon persistence

## Verdict
`AUTHENTIC_DUAL_KAPPA_PARAMETER_SENSITIVITY_AND_EXTENDED_MODEL_CLOCK_IVPS__SAMPLED_GEOMETRY_PERSISTENT`

## Completed integrations
26

## Provenance
v18.51s: SHA256_AUTHENTICATED; native v18.52l: SHA256_AUTHENTICATED

## Scientific interpretation
The 12.8–25.6 interval is independently checked against the authenticated v18.53h fix1 and original native trajectory.
The 25.6–51.2 extension, if computed, is a continued ORIGINAL-SOURCE forcing protocol, NOT independently authenticated historical observations.
A source-matched IVP may diverge parametrically without invalidating the original model. Report angles, derivative conditioning and time transport independently.
Principal-plane projectors are meaningful away from singular-value degeneracy; Schur bases themselves have arbitrary gauge.
Time-integrated bivector compensation does not certify instantaneous cancellation, global stability, or hinge Noether charge.

## Summary
{
  "number_parameter_values": 11,
  "number_seeds": 2,
  "number_integrations": 26,
  "number_parameter_plane_comparisons": 120,
  "number_horizon_windows": 132,
  "number_projector_derivatives": 48,
  "number_time_transport_comparisons": 12,
  "weak_geometry_wider_kappa_cases": 0,
  "max_charge_drift": 1.021405182655144e-14,
  "max_authentic_rhs_error": 2.0660556598883772e-15,
  "max_numeric_refinement_error": 1.1401990462900358e-13,
  "max_historical_shared_replay_error": 9.059419880941277e-14,
  "max_historical_53h_residual_replay_error": 9.992007221626409e-16,
  "max_quadrature_defect": 0.00011407152297471188,
  "min_principal_gap": 0.8274307408896822,
  "source_history_expires_xi": 25.6,
  "continued_clock_end_xi": 51.2
}
