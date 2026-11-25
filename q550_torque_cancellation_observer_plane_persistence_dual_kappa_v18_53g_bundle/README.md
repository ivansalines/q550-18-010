# Q550 — v18.53g: Torque-Cancellation Observer-Plane Persistence & Dual-κ₀ Audit

This reproducible **diagnostic** examines whether the radial-potential versus phase-coupling torque cancellation from authenticated v18.53f persists in seven predefined observer planes, for both historical κ₀ prescriptions and their shared midpoint.

**Key distinction:** the source v18.52l trajectories were integrated with **one shared κ₀**. Evaluating the authenticated RHS at κ_autonomous and κ_original on these *same archived trajectories and same historical external phases* gives a **frozen-state counterfactual vector-field comparison**, not new self-consistent dynamical trajectories. No new IVPs are performed.

## Exact provenance and operational requirements

- `reference/v18_53f_authenticated_report.json`, `reference/v18_53f_authenticated_native_summary.json`, and the previous output manifest are retained with SHA-256 hashes and cross-checks.
- The historical v18.51s nested parent source ZIP is required for source-level torque attribution. It is authenticated against the pinned SHA-256 before its original generator is staged. Absent parents yield a diagnostic sendback but **no false certificate**.
- The historical v18.52l 10D native NPZ file must match its exact SHA-256. It is located automatically in the same server folder hierarchy (or passed explicitly).
- Source comments, protocol, and reports are in English.

## Commands

```bash
cd /root/elab1/hysteresis
unzip q550_torque_cancellation_observer_plane_persistence_dual_kappa_v18_53g_bundle.zip
cd q550_torque_cancellation_observer_plane_persistence_dual_kappa_v18_53g_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent --require-native 2>&1 | tee ../v18_53g_server.log
```

Manual source paths: `--parent-v1851s /path/to/..._server_sendback.zip --native-npz /path/to/extended_horizon_10d_states.npz`.

Outputs: 420-row projection×κ₀×series table, invariant integrated bivector table, signed-time/local opposition diagnostics, baseline anchoring, RHS spot checks, plots, scientific report and full sendback archive.

## Classification policy

- A cancellation score over 0.9 is *defined* as strong cancellation in the two signed integrated projected channels; a score is **not** proof of a time-local or invariant physical cancellation.
- The coordinate-free integrated bivector comparison provides a separate observer-independent condition.
- A `PERSISTENT_ACROSS_PREDECLARED_PROBES` label requires thresholds to pass for **every** predefined plane and both κ prescriptions, not only the original plane. Otherwise report partial persistence and failing angles explicitly.
- A threshold pass is numerical and bounded by the explored source states/clock; it is not a global law.
- There is no hinge Noether-charge or autonomous-rule trajectory certificate here.
