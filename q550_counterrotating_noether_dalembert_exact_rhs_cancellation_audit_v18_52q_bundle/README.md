# Q550 — v18.52q: Counter-Rotating Noether–d’Alembert Exact RHS Cancellation Audit

**Objective:** Verify from the authentic nonlinear Q550 angular equations whether
`C = exp(2*l1)*omega1 - exp(2*l2)*omega2` is an **exact first integral**, not a
fortuitous numerical invariant of the observed orbit. Distinguish dynamical
conservation, group equivariance, a Noether–d’Alembert interpretation, a
conditional noncollapse barrier, and topological winding (not established).

## Historical provenance

- The original 2026-08-31 `test_support_q550_oriented_coupled_bulge_cone_audit_v1.py` is included **unchanged**. It carries the explicit 11D nonlinear RHS and support-clock law.
- The independent 2026-08-11 `test_support_q550_support_ode_interval_residual_pilot_v1.py` is included unchanged and confirms the same angular structure.
- Three authentic v18.52g/v18.51t loader and kinetic-embedding source files are preserved, for loading the exact `structural_jet` engine from the v18.51s server sendback.
- The v18.52l 10D source is SHA-256 pinned. It is referenced **in place**, not duplicated inside this lightweight ZIP.
- The v18.51s sendback archive is SHA-256 pinned. The server audit will not silently substitute a different generator.

## Run

```bash
cd /root/elab1/hysteresis
unzip q550_counterrotating_noether_dalembert_exact_rhs_cancellation_audit_v18_52q_bundle.zip
cd q550_counterrotating_noether_dalembert_exact_rhs_cancellation_audit_v18_52q_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_52q_server.log
```

`run_on_server.sh` verifies source SHA-256, runs independent tests, evaluates the
archived nonlinear field, evaluates all 20 native traces, and stages the **authentic
v18.51s→v18.51r→... structural_jet** for independent off-trajectory symmetry and
conservation probes, without re-integrating the Q550.

If the v18.51s server sendback has been moved, pass a full path:
`bash run_on_server.sh --source-s-path /absolute/path/to/...v18_51s_bundle_server_sendback.zip`

Development-only: `python3.13 scripts/run_audit.py --mode archived` skips the
v18.51s live-engine probe, but its report MUST NOT be presented as full-engine
verified. All generated scripts, reports, and source reference files are included
in `...v18_52q_bundle_server_sendback.zip`.

## Decision boundary

**Exact algebraic cancellation** is a theorem *for the documented RHS*. The
**actual engine comparison** verifies that this is the model used in the present
experiment. `p_theta1-p_theta2` is a generalized Noether–d’Alembert charge when
the angular forces are the Euler–Lagrange forces derived from the stated kinetic
and interaction terms, with remaining nonconservative forces acting in radial
and auxiliary coordinates. It is **not a quantized topological charge**.

A nonzero conserved `C` and an independently guaranteed upper kinetic-energy
bound forbid **simultaneous collapse of both absolute radii**. They do NOT, by
themselves, forbid individual-radius collapse, relative separation collapse,
or future unbounded kinetic-energy growth.
