# Q550 — v18.53s
## Quadrature-Order, Cumulative-Primitive Bias & Same-Grid Physical-Step Separation Audit

**Status:** Server distribution. The full physical audit has **not** been executed locally.

This version addresses the numerical difference left open by v18.53r: 4x and 2x sources at common grid nodes differ at ~1e-13, while cumulative primitives and valley geometry retain larger, systematically converging differences. It separates *quadrature sampling* from *physical 11D integration step/tolerance* without altering the underlying Q550 equations.

### Inputs and full genealogy

- `inputs/v18_53r_authenticated_server_sendback.zip` is an immutable pinned physical checkpoint (12 new IVPs); it contains the entire nested q → p → o → n/m/l/k/j/h model ancestry and original executable scripts.
- The bundle does not duplicate the two large native historical input sources. Full physical execution finds the user's original hashed files under `/root/elab1/hysteresis` (or uses the supported explicit CLI overrides).
- All scripts, scientific criteria, source hashes, tests, manifest checks, and regenerated outputs live together in the archive.

### Independent numerical interventions

1. **Quadrature only, one held-fixed 4x physical trajectory:** composite trapezoid on nested 1x/2x/4x nodes; composite Simpson on 2x/4x nodes; Richardson-style contrasts as *indicators*, not accuracy guarantees; exact source-PCHIP integral at the same nodes.
2. **Telescoping split:** `T4(g4)-T2(g2) = [T4(g4)-T2(g4 restricted)] + [T2(g4 restricted)-T2(g2)]`. This exactly accounts for integration density and the change of source values at common nodes. It is not a causal physical decomposition.
3. **Same 4x grid, tighter physical integration:** 6 new Q550 11D DOP853 replays with the same `t_eval` samples as v18.53r, a halved `max_step`, and 10x tighter rtol/atol. This independently probes changes induced by numerical 11D evolution while holding the quadrature rule and output grid fixed.
4. **Oriented valley diagnostics:** compare the archived cumulative primitive with reconstructed trapezoid and source-integrated PCHIP primitives, plus the new tighter 11D replay. At the same origins and widths, report the minimum 10D residual, the width path, direction speed, direct angular arclength, and integrated oriented 45D area.

### Execute on the user's server

```bash
cd /root/elab1/hysteresis
B=q550_quadrature_order_cumulative_primitive_bias_physical_step_separation_v18_53s_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
 bash "$B/scripts/run_on_server.sh" \
 --require-parent --require-native --require-dynamic \
 2>&1 | tee v18_53s_server.log
```

Fast tests and authenticity checks without new IVPs:

```bash
python3.13 "$B/scripts/run_audit.py" --preflight
python3.13 "$B/scripts/run_audit.py" --smoke
```

Optional diagnostic-only execution (does **not** certify step separation):

```bash
PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" --quadrature-only
```

If using the diagnostic-only path, archive the resulting sendback under an explicitly separate name: the full server execution must not be confused with postprocessing.

### Outputs after physical execution

- `results/quadrature_order_bias.csv`
- `results/physical_step_reference.csv`
- `results/same_grid_physical_source_comparison.csv`
- `results/valley_geometry_quadrature_bias.csv`
- `results/oriented_valley_quadrature_comparison.npz`
- `results/quadrature_geometry_separation_*.png`
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`
- A generated, independently CRC and SHA-verified `...v18_53s_bundle_server_sendback.zip`.

### Scientific boundaries

All comparisons use the Euclidean 10D bivector convention of v18.53r. The formal clock continuation beyond ξ=25.6 is retained. Neither empirical quadrature order, same-grid tighter ODE replay, minimum of the 10D residual nor oriented area establishes a true full-state return, material spin or material holonomy.
