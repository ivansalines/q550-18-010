# Q550 — v18.52s: Energy–Work–Dissipation Balance & Global Noncollapse Bound Audit

**Do not run another long IVP.** This is a reconstruction from the archived exact RHS, and the authenticated v18.52l 10D trajectories.

From `/root/elab1/hysteresis`:

```bash
unzip q550_energy_work_dissipation_global_noncollapse_bound_v18_52s_bundle.zip
cd q550_energy_work_dissipation_global_noncollapse_bound_v18_52s_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52s_server.log
```

Expected native 10D input (not embedded, SHA-256 checked):
`../q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz`.

The authenticated global phase is NOT inside the 10D NPZ: it is reconstructed from `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`, SHA-256 `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`, using original `ctx[theta](event+xi)` and the v18.51t genealogy. The 11th native state column is material tau, NOT phase Theta. Without the parent archive the energy ledger is not run and the report is explicitly partial. Use `--parent-v1851s PATH` if necessary.

If relocated, pass `--states /absolute/path/to/extended_horizon_10d_states.npz`; if missing, symbolic/off-orbit still run and report is explicitly partial. Use `--require-traces` for fail-closed behavior.

Outputs: `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `results/native_20_series_energy_ledger.csv`, optional PNGs and symbolic result. `scripts/package_sendback.py` packages outputs and all scripts.

**Scientific scope:** Exact energy balance, a rigorous *finite-physical-time* Gronwall bound excluding finite-time collapse of absolute radii when e_reserve>0, and 20-series numerical ledger consistency. This **does not** prove a uniform-in-time bound or noncollapse of the relative 5D perturbation vector.
