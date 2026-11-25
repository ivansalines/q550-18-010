# Q550 — v18.53t
## On-the-Fly 21D Primitive Integration, Adaptive-Step Coupling & Full-10D Valley Audit

**Status:** server-ready distribution; physical experiment is **not** executed in this distribution.

This test changes the **numerical integration protocol**, not the original Q550 physical equations. It appends ten passive source integrals to the original **11D Q550** system and integrates the augmented 21D system with DOP853. The 10D source is the historical Euclidean bivector **radial potential + phase coupling**, expressed per unit of xi clock (with the original material-time conversion). The accumulator is mathematically passive: it does not appear in the physical 11D RHS.

## Scientific question

The v18.53s audit showed that changes in sampled quadrature (~2.33e-6) dominated differences from tightening the 11D physical timestep (~2.67e-13) at the tested resolution. Does direct on-the-fly integration of the 10D primitive preserve the positive near-return minimum and oriented 10D transport? How much does adding an accumulator change adaptive time steps and, therefore, the physical numerical trajectory?

## Authenticated genealogy

`inputs/v18_53s_authenticated_server_sendback.zip` is a SHA256-pinned original server sendback. It embeds the r→q→p→o→n→m→l→k→j→h sources, data and manifest genealogy. The original source is not modified. Historical v18.51s physical parent and native NPZ are large external model checkpoints; they must be available on the user's server and pass their original SHA256 checks. No synthetic input is substituted.

## Three integration configurations

1. **11D physical-only DOP853** at fixed physical tolerances, `max_step=0.004375` in xi. Save its accepted-step mesh and dense output.
2. **21D fully augmented DOP853** with the identical 11D RHS and ten passive integral coordinates. Save the accepted-step mesh, 11D states, accumulated 10D primitive, source and dense output. The additional coordinates may affect acceptance via the solver's composite error norm: this influence is explicitly measured, never assumed absent.
3. **10D accumulator on the *locked* 11D dense trajectory**, integrated independently with DOP853 without feedback on the physical 11D solver. This is a control for accumulator error on one already-computed numerical trajectory.

The on-the-fly 21D primitive is also compared with the original s sampled-source PCHIP integral and a same-grid 4x trapezoid. All reconstructed valley geometries reuse the pre-registered original origin/width domain and the same Euclidean 10D norm.

## Reproducible server command

```bash
cd /root/elab1/hysteresis
B=q550_on_the_fly_21d_primitive_integration_adaptive_step_coupling_full_10d_valley_v18_53t_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
  bash "$B/scripts/run_on_server.sh" \
   --require-parent --require-native --require-dynamic \
   2>&1 | tee v18_53t_server.log
```

The script automatically builds and verifies `q550_on_the_fly_21d_primitive_integration_adaptive_step_coupling_full_10d_valley_v18_53t_bundle_server_sendback.zip` only after all 18 required integrations have completed. Failure is explicitly recorded in `results/report.json` and cannot be packaged as success.

### Preflight without physical IVPs

```bash
PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" --preflight
```

The preflight authenticates s/r and their ancestry, checks source manifests and executes only synthetic and source-construction tests. Do not interpret its success as a physical result.

### CLI overrides for original native checkpoints

If historical checkpoints are not at the standard server locations, pass `--parent-v1851s /absolute/path/to/original_v1851s.zip --native-npz /absolute/path/to/original.npz`. Both hashes are checked against the model history.

## Outputs

- `results/physical_21d_vs_11d.csv` — state/source/primitive comparisons and accepted-step mesh diagnostics.
- `results/valley_geometry_21d.csv` — full-10D minimum, center directional speed, direct angular length and 45D-oriented polygon.
- `results/augmented21_primitives_and_oriented_valleys.npz` — 21D physical states, primitives, 10D sources, oriented direction samples and direct dense-output checks.
- `results/augmented21_valley_original.png` and `results/augmented21_valley_autonomous.png` — graphics.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json` — scientific decision and provenance.

## Scope and limitations

The physical model is not redefined; the 21D system merely appends the integral of the previously defined 10D source. The solver's adaptive step decisions need not remain identical. Same-grid dense output is a numerical interpolant, not an exact solution. The oriented 10D bivector residual is **not** a full material-state recurrence or a material holonomy. The studied late endpoints are a **formal clock continuation beyond xi = 25.6**. No new physical invariant or global minimum is presumed.
