# Q550 — v18.52p
## Canonical Hinge, Noether-Torque & Conditional Noncollapse Barrier Audit

This is a **falsifiable postprocessing audit**, NOT an a priori conserved-charge declaration. It separates (a) exact canonical coordinates, (b) physical absolute angular-momentum *kinetic proxies* derived from the exact Q550 kinetic embedding, (c) the relative separation angular momentum between two distinct nonlinear trajectories, (d) torque and radial identities, (e) conditional barriers to planar/5D collapse, and (f) open-path angular winding.

**Critical distinction:** a symplectic change of chart preserves the expression `J=Q1 P2-Q2 P1` identically, but only a continuous symmetry of the *actual model action* implies `dJ/dtau=0` by Noether. The archived Q550 field includes a driven nonlinear 10D flow and a material clock; a complete Hamiltonian/action with conjugate momenta has NOT been reconstructed here. No conservation or topological claim is assumed from its projected geometry.

## Sources and trust

These sources MUST already exist in `/root/elab1/hysteresis`:

```
q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/
  results/extended_horizon_5d_vectors_and_traces.npz
  results/extended_horizon_10d_states.npz
```

Both authentic v18.52l SHA-256 values are hard-pinned, and the historical v18.49h `kinetic_metric_v18_49h.py` code is preserved, with its SHA-256. `reference/PARENT_v18_52o_report.json` anchors the previous authenticated full-path 3:1 and hinge observer geometry. **No nonlinear IVP is rerun**.

The original four observed Cartesian point coordinates are still not provided. The user-observed slope `m=-0.308` is therefore a *pinned observational hypothesis*, not a regression fitted from those points or a discovered dynamic invariant. `atan(m)` is approximately **−17.119°**, not exactly −17.04°. The ±17.04° paper tilt around the axis is an independent observation construction.

## Run on the user server

```bash
cd /root/elab1/hysteresis
unzip q550_canonical_hinge_noether_torque_and_noncollapse_barrier_audit_v18_52p_bundle.zip
cd q550_canonical_hinge_noether_torque_and_noncollapse_barrier_audit_v18_52p_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh 2>&1 | tee ../v18_52p_server.log
```

Workflow: verify ancestry → synthetic falsification tests → genuine SHA-verified preflight → **full physical 5D and native 10D** analysis → output hashes → complete ZIP sendback.

A development-only `--sample-stride 128` and `--allow-missing-10d` are available for local pipeline tests, but never produce scientific-certified/full-resolution reports. The server defaults to `--sample-stride 1` and REQUIRES both authentic sources.

If needed, override source locations by passing arguments to `bash run_on_server.sh --source5 /absolute/path.npz --source10 /absolute/path.npz`. Adjust `--output-stride 32` to decimate only the exported CSV, NOT metric calculations.

## Main output files

- `results/report.json` — all source hashes, branchwise 5D and native 10D summary, 2D/5D angular decomposition, angular torque and collapse status.
- `results/SCIENTIFIC_REPORT_SERVER.md` — evidence tables and candid interpretation.
- `results/relative_5d_hinge_kinematics_decimated.csv` — explicit samples of torque and radial competition.
- `results/observer_gate_kinetic_torque_events.csv` — ±17.04° observer sections tagged with angular and radial data.
- `results/native_absolute_vs_relative_J_matched_time.csv` — **direct check that `J_kick-J_reference` is not `deltaQ cross deltaV`**.
- `results/relative_hinge_torque_and_collapse_balance.png` — relative material hinge J, torque, radius and acceleration balance.
- `results/native_absolute_kinetic_moment_comparison.png` — native absolute first-rotor kinetic J curves and their difference.
- `OUTPUT_MANIFEST.json` and `SOURCE_MANIFEST.json`, `reference/` and `scripts/` — runnable genealogy.

The sendback archive is created in the parent folder with a matching full descriptive root:

`q550_canonical_hinge_noether_torque_and_noncollapse_barrier_audit_v18_52p_bundle_server_sendback.zip`.

## Scientific interpretation limits

- `Jrelative=q1*v2-q2*v1` uses differences between *two separate experiments* and is **not** a derived conserved Noether generator of their physical evolution.
- `Jabsolute=r1²*w1` is the kinetic momentum conjugate to the first rotor angle for the `T=1/2||V||²` kinetic term. Full canonically conjugate momentum and original force symmetry require the actual action/force terms; numerical angular balance **may falsify conservation**, but does not prove Noether symmetry.
- `J' = J` under the same orthogonal planar rotation of `(Q,P)` is an exact algebraic result even when `J(tau)` is not constant.
- `T_plane = 0.5*rdot² + J²/(2r²)` and `rddot = J²/r³ + (q·a)/r` are exact kinematic identities wherever `r>0`. A **dynamical centrifugal barrier against collapse** requires a nonzero conserved charge, appropriate kinetic metric, and an energy or velocity bound for the future solution. A positive sampled radius is an observation, not a global theorem.
- Unwrapped angular turns for the open time interval are not topological winding integers.
