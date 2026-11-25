# Q550 — v18.53r: Denser Physical-Source Sampling, Full-10D Valley & Oriented Transport Convergence

This audit is the first new physical Q550 replay after v18.53o in this sequence. It checks whether the source-constrained valley geometry seen in v18.53q is stable when the **original 11D model** is integrated again on 2× and 4× denser sample grids with reduced adaptive DOP853 maximum steps. The six original seed/κ₀-rule combinations are retained exactly. **All computations of the full trajectories run on the user's server, not locally.**

## Scientific measurements

- Authenticate the entire v18.53q → v18.53p → v18.53o and physical ancestor chain. Preserve it byte-for-byte in `inputs/`.
- Retrieve the historical physical model implementation by loading the authenticated v18.53o genealogy; never substitute a new ODE or fabricate Q550 data.
- Reintegrate the 11D physical state over exactly the old ξ interval, sampling every historical grid interval at 2× and 4× density, and reducing the DOP853 `max_step` proportionally. Tight tolerances are selected from the authenticated parent. Increasing `t_eval` alone would not certify denser time stepping.
- Extract the actual force-derived radial/coupling bivector sources, physical clock-rate factor, and cumulative primitive on each new grid. Record absolute and relative differences of C and of the source at the original knots.
- Reconstruct 10D unit direction, width-optimized valley, direct geodesic length, chain-rule length, and the 45-component oriented wedge. Compare Hermite source-constrained and PCHIP on the 1×, 2×, 4× genuine model grids. Compare 1× baseline using archived original physical samples (not a new integration).
- Report per-case convergence and a conservative scientific decision. No full-state material recurrence, physical spin, or holonomy is inferred from bivector coordinates.

## Server run

```bash
cd /root/elab1/hysteresis
B=q550_denser_physical_source_sampling_full_10d_valley_oriented_transport_convergence_v18_53r_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 bash "$B/scripts/run_on_server.sh" --require-parent --require-native --require-dynamic 2>&1 | tee v18_53r_server.log
```

The original external 11D parent and native Q550 source files remain required on the server and are checked against their historically pinned SHA-256 values. CLI overrides `--parent-v1851s PATH` and `--native-npz PATH` are available. A local, no-IVP audit is `python3.13 scripts/run_audit.py --preflight`. `python3.13 scripts/run_audit.py --smoke` reuses the parent traces only.

The output is `q550_denser_physical_source_sampling_full_10d_valley_oriented_transport_convergence_v18_53r_bundle_server_sendback.zip`. It contains exact Python sources, results CSV/NPZ, plots, report, ancestor archives, manifests and provenance. All comments, documentation and figures are in English.

## Limits

These are denser **numerical solutions of the same historical 11D dynamical equations**, not experimentally independent physical trajectories. The force/source histories past ξ=25.6 still belong to the **formal clock continuation** of the model. The Euclidean norm in 10D bivector coordinates is not a physical metric on the complete state. The derivative at a knot and the width optimizer may still introduce numerical sensitivities; report them rather than hiding them.
