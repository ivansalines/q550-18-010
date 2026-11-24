# Q550 — v18.52v-r1: Forward–Return Signed Drive, Clock Orientation, Fixed-Kappa & Energy Quadrature Audit

Corrected reproducible companion to v18.52v. Does NOT reconstruct an unobserved forward orbit.

## Corrections justified by authenticated ancestor
1. `cubic_phase_excursion` differentiates to **minus** `cubic_clock_sign` on the regular historical ramp/return intervals. The original `theta(t)` itself includes a minus; therefore `dtheta_return/dxi = +omega*cubic_clock_sign`, exactly the signed original RHS's eleventh component. The v18.52v diagnostic compared the opposite signs.
2. The v18.52l native long-horizon 10D trajectories were evolved at the common, fixed `kappa_native=(kappa_original+kappa_autonomous)/2 = 0.1666331385421536`. The original v18.52v incorrectly evaluated each native seed with its historical seed kappa, creating a model-parameter mismatch. This version FAILS CLOSED unless the authenticated v18.51s fixed-kappa provenance matches.
3. Evaluate **all 20 native series**, including refined, at strides `1 2 4 8 16 32`. Full res means 65,537 samples per series. The phase rate is derived from the authenticated source, not chosen by minimizing the closure error; the numerical derivative is an independent check.
4. Compare an intentionally wrong-kappa on-path energy ledger for falsification, not parameter fitting.

## Server run
```bash
cd /root/elab1/hysteresis
unzip q550_forward_return_signed_drive_orientation_kappa_reconciliation_and_energy_quadrature_v18_52v_r1_bundle.zip
cd q550_forward_return_signed_drive_orientation_kappa_reconciliation_and_energy_quadrature_v18_52v_r1_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent --require-native 2>&1 | tee ../v18_52v_r1_server.log
```

Inputs from server (not duplicated):
- `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` — SHA-256 in `inputs_config.json`.
- `q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz` — SHA-256 in config.

`results/` contains server scientific report, signed ledger CSV for every stride, kappa falsifier CSV, authenticated source extracts, source-protocol CSV, figure PNGs and SHA-256 output manifest. The server packaging script produces a matching-root `*_server_sendback.zip` even if a required ancestor is missing (nonzero exit code in strict mode).

No new native Q550 integrations. Historical `make_context` may stage its inherited short ancestor IVP as in v18.52v.
