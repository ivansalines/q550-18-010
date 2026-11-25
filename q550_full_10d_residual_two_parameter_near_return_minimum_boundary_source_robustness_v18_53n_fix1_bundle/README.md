# Q550 — v18.53n-fix1 Full-10D Residual Two-Parameter Near-Return Minimum & Boundary-Source Robustness Audit

## Scientific question

Does the observer-independent **full 10D bivector norm** of the telescoped three-window residual have a stable two-parameter minimum in the region identified by v18.53m, or was apparent closure limited to fixed 2D SVD projections?

This audit deliberately does **not** equate a small integrated bivector residual with full 11D phase-space return, physical holonomy, dynamical recurrence, or material memory.

## Inputs and provenance

- Embeds the **byte-identical authenticated v18.53m-fix1 server sendback**, including v18.53l and its nested genealogy, Python scripts and generated results.
- Requires the same two original server physical source inputs: ancestral v18.51s parent source and native Q550 NPZ. SHA-256 values are pinned in `inputs_config.json`; absent or wrong inputs cause a hard failure.
- Parent's exact 36 projected roots are verified against the newly integrated full 10D residual before optimizer results may be accepted.
- Script comments and report text are English. No heavy physical Q550 integration occurs during local preflight.

## Server run

```bash
cd /root/elab1/hysteresis
B=q550_full_10d_residual_two_parameter_near_return_minimum_boundary_source_robustness_v18_53n_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
 bash "$B/scripts/run_on_server.sh" \
 --require-parent --require-native --require-dynamic \
 2>&1 | tee v18_53n_server.log
```

The script creates `$B`_server_sendback.zip next to the extracted root.

To verify genealogy without integrating physics:

```bash
python3.13 "$B/scripts/verify_manifest.py"
python3.13 "$B/tests/test_synthetic.py"
python3.13 "$B/scripts/run_audit.py" --preflight
```

## Experiments

- Two physical seeds, 3 κ₀ rules: 6 primary 11D integrations + 4 tighter integrations.
- Full 10D Frobenius-isometric net norm, **no projection chosen to minimize**.
- Grid: 53 origins × 43 widths × 6 physical cases = **13,674** sampled objective values.
- Bounded multistart 10D least-squares minimizations, KKT gradients, 10D Jacobian singular values and conditioned local Hessian of the interpolant.
- Compare against v18.53m's 36 projected roots, calculate full net gain.
- Derive boundary source Jacobian `J=[g_end-g_start, 3 g_end]` and compare with physical source interpolation.
- Test primary vs every-other-sample cumulative quadrature, PCHIP interpolation sensitivity, and tighter physical IVPs.
- Report boundary minimizers, trace knot clearance, and local stationarity, and avoid global-optimum claims.

## Outputs

`results/` will contain the 13,674-row grid, all multistart solutions, deduplicated local minima, six best minima, parent-root replay, quadrature/tight-IVP sensitivity, authenticated physical run metadata, plots and English scientific report. `OUTPUT_MANIFEST.json` hashes the results.

## Conditions of scientific interpretation

- Three interior window boundaries cancel exactly in the sum. **Only two time endpoints determine the net bivector**.
- An exact projected root in any 2D plane does not imply a full 10D zero.
- Even a full 10D zero would not establish material return of all physical dynamical variables.
- A locally quadratic objective between sampled knots is a property of **linear source interpolation**, not a proof of smooth physical second variation.
- The clock past ξ=25.6 is the formal continuation of the archived model, not an independently measured trajectory.
- The protocol explicitly allows a null result and documents non-convergence.

## v18.53n-fix1 notes

The authenticated original v18.53n server sendback is preserved unchanged under
`provenance/`. It stopped after one physical integration with
`NameError: name 'net_j' is not defined` in `endpoint_source_diff`.

Fix1 changes only `net_j` to the already defined `net_jac` in that helper.
It adds a synthetic test that exercises the endpoint/source comparison and
keeps all Q550 equations, parameter bounds, stopping rules, and physical
source inputs unchanged. All physical results still require server execution.

The root directory and ZIP filename end in `v18_53n_fix1_bundle`. The
returned server archive ends in `v18_53n_fix1_bundle_server_sendback.zip`.
