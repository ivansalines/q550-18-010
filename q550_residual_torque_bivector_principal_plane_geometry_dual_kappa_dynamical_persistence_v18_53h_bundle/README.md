# Q550 — v18.53h: Residual-Torque Bivector Principal-Plane Geometry & Dual-κ₀ Dynamical Persistence Audit

This executable genealogy follows the **entire physical 5D integrated residual torque bivector** of authenticated v18.53g and performs a genuinely **independent original/autonomous κ₀ integration**, rather than just modifying κ in the RHS while holding the path fixed.

## Primary research questions

1. What are the two principal orthogonal planes, singular intensities, and kernel of the 5D skew residual `M_radial + M_phase`? Do they remain identifiable at different times? (Near-degeneracy is explicitly flagged.)
2. How much of the 5D cancellation is invariant, and how much arises only in the original observer plane?
3. Starting at the **same authenticated state at ξ = 12.8**, do separate integrations for original and autonomous κ₀ exhibit comparable residual-bivector geometry through ξ = 25.6?
4. Can a shared-κ₀ trajectory be reproduced from the original native initial state, and are independently integrated results converged as tolerances tighten?

**Important historical limit.** The separate κ runs are self-consistent *over the recalculated late-window*, with identical initial conditions and authenticated external clock. The pre-ξ=12.8 historical ramp is not recalculated under the two prescriptions. No full historical bifurcation, whole-orbit conservation of hinge moment, or unique canonical hinge symmetry is certified.

## Standard server command

```bash
cd /root/elab1/hysteresis
unzip q550_residual_torque_bivector_principal_plane_geometry_dual_kappa_dynamical_persistence_v18_53h_bundle.zip
cd q550_residual_torque_bivector_principal_plane_geometry_dual_kappa_dynamical_persistence_v18_53h_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent --require-native --require-dynamic 2>&1 | tee ../v18_53h_server.log
```

Manual input overrides (if automatic discovery fails):

```bash
PYTHON_BIN=python3.13 bash run_on_server.sh \
 --parent-v1851s /root/elab1/hysteresis/q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip \
 --native-npz /root/elab1/hysteresis/q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz \
 --require-parent --require-native --require-dynamic
```

A 97 MB native NPZ and the historical v18.51s nested genealogy are external **content-hash authenticated inputs**; they are intentionally not duplicated inside this small portable ZIP. Their SHA-256 values are pinned in `inputs_config.json`. The local smoke run without those inputs produces a marked *partial diagnostic sendback* and nonzero return code under the `--require-...` flags.

## Outputs

- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`: scientific classification; provenance and validation failures if any.
- `results/dynamically_integrated_dual_kappa_orbits.csv`: actual independent κ dynamical integrations, source agreement, convergence, conservation, residual bivector metrics.
- `results/bivector_principal_window_geometry.csv`: per-seed/per-κ complete and four window principal singular-value diagnostics.
- `results/paired_kappa_principal_planes.json`: angles between principal planes and residual bivector comparisons at fixed seed and fixed ξ endpoints.
- `results/paired_kappa_dynamical_geometry.csv`, `results/window_cancellation_persistence.csv`: parameter differences, invariant persistence across intervals.
- `results/baseline_53g_reconciliation.json`: 53g integrated-bivector replay using original native arrays before any new orbit integration.
- `results/dual_kappa_dynamic_radius_and_charge.png`: physical radius and J conservation of the separately integrated late-window orbits.
- automatic fully hashed `..._server_sendback.zip` in the parent directory.

The old v18.53g **same-state** scores (~79.5% invariant, >96% in an original plane) are **inputs to the question, not expected answers** for the new dynamical paths. The script never overwrites them as if they were new dynamic results.
