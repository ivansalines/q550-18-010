# Q550 — v18.52g: Inward–Kinematic Competition & Material-Time Horizon Audit

**Scope:** genuinely extend the historical nonlinear 10D Q550 finite-kick trajectories beyond the v18.52f/v18.52e event-clock horizon `xi<=0.2` and find whether the inward relative-acceleration projection can overcome the positive kinetic term. Three distinct events are tested at common *material time*: the convexity zero `D''=0`, the separation-rate zero `D'=0`, and a maximum of squared position separation (positive-to-negative `D'` crossing).

The server-run audit **does not extrapolate the old f plot**. It re-stages the original authenticated v18.51s→r→q→19-ancestor 10D vector field, recreates the exact five-dimensional Cartesian kinetic embedding from v18.49h, reintegrates baseline and lateral finite-kick branches with **unchanged** midpoint κ and physical parameters, and aligns branch observations by a numerically integrated material clock. It compares the newly calculated first 0.2 of the trajectory with the archived history before interpreting the extension.

## One-time dependency: original authenticated v18.51s parent

The source bundle includes the complete exact **v18.52f server sendback** (which embeds the v18.52e and v18.52d sendbacks), reference scripts, and historical physical-state arrays via that chain. The v18.51s archive needed to restore the full original 10D engine was **not contained in the user-provided f sendback** and is therefore a required external source already used for earlier server executions:

`q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`

Expected SHA-256: `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`.

Place it at `/root/elab1/hysteresis/` (the standard server workspace), next to the extracted bundle, or pass its absolute path using `--source-s-path`. If it cannot be found or its hash differs, the server run **stops**; no unverified synthetic model is substituted.

## Server execution

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_inward_kinematic_competition_material_time_horizon_v18_52g_bundle.zip
cd q550_inward_kinematic_competition_material_time_horizon_v18_52g_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  --xi-end 0.8 --samples 2049 --refine \
  2>&1 | tee ../v18_52g_server.log
```

If the original field becomes ill-conditioned beyond the historical domain, **do not** change the RHS or silently extrapolate the data. Try a shorter physical horizon, e.g. `--xi-end 0.4`, and retain the error log. For a later, independent exploratory horizon use `--xi-end 1.2` only after the first extension is checked. `--max-step` defaults to `0.00125`; the `--refine` run uses `0.000625` and a fourfold tightening of relative and absolute tolerances. Refinement is enabled by default; `--no-refine` is for a preliminary run only.

A fast metadata and exact-parent-only local test (no nonlinear model integration) is also available:

```bash
python3.13 tests/test_synthetic.py
python3.13 tests/test_mock_clock_ivp.py
python3.13 scripts/run_audit.py --mode preflight
```

## Outputs

- `results/report.json`: signed physical diagnostics, source lineage, per-kick status, event brackets, convergence checks, all caveats.
- `results/material_horizon_ledger.csv`: per-section physical material time, state difference, inner/outer competition and drive-phase mismatch.
- `results/material_horizon_events.csv`: located C and G roots and D maxima, including an empty-header file when none are detected.
- `results/kick_horizon_summary.csv`: per-kick condensed summary.
- `results/material_horizon_traces.npz`: reproducible physical scalar time series.
- `results/extended_10d_nonlinear_states.npz`: full ten-dimensional seed and kicked physical states plus material clock, for the requested event-clock grid (coarse and refined when enabled).
- `results/inward_kinematic_material_horizon.png`: diagnostic figure, with historical horizon marked.
- `results/SCIENTIFIC_REPORT_SERVER.md`: report without substituting successful execution for a positive physical result.
- `OUTPUT_MANIFEST.json`: SHA-256 of generated files.
- `../q550_inward_kinematic_competition_material_time_horizon_v18_52g_bundle_server_sendback.zip`: exact packed sendback including executable generator and input provenance.

Python >=3.10 with NumPy, SciPy, Matplotlib. The scripts and document comments are **in English**. All expensive nonlinear integrations must run on the user's server.

## Scientific stopping rules

A zero of `D''` is **not** a position return. A zero of `D'` is not a maximum unless it crosses from positive to negative. The material time horizon ends at the earlier of the reference or kicked trajectory endpoints, so neither is extrapolated. Even a numerically verified maximum remains a finite-horizon result in the **specified kinetic Cartesian metric**, not a proof of asymptotic recovery, conserved charge, or gauge invariance. Dense-output numerical roots are not rigorous interval certificates and may miss tangencies or roots occurring between grid checks.
