# Q550 v18.53z — 11D Tangent Generator, Complement-4D Rotation, Lift Sensitivity & Gauge-Covariant Transport

This audit starts from the immutable, authenticated **v18.53y server sendback** and its complete ancestry through x/w/v/u/t. It asks whether a **local, finite-lag linearization of the authenticated Q550 11D equations** selects a complementary rotation beyond the minimum-turn kinematic convention in y.

**Important:** The model determines an 11D tangent generator, but projecting it to a 5D configuration rotation is **not uniquely determined**. The 5D-to-11D horizontal lift and finite-lag horizon are explicit. A nonzero complement does not by itself demonstrate a physical material connection, holonomy, or 11D recurrence.

## Server execution

```bash
cd /root/elab1/hysteresis
B=q550_11d_tangent_generator_complement_4d_rotation_lift_sensitivity_gauge_covariant_transport_v18_53z_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
bash "$B/scripts/run_on_server.sh" --require-parent 2>&1 | tee v18_53z_server.log
```

`--preflight-only` authenticates the sources and parents and runs synthetic tests. `--smoke-only` executes a genuine historical case at three lambda values without writing scientific reports. The default `--require-parent` executes the full audit and packages the sendback.

## Computation

- No new **nonlinear** physical 11D/21D integrations.
- Six authenticated historic physical trajectories, two epochs, 17 values of metric lambda, three temporal window assignments.
- 204 parent paths × 3 start/mid/end lifts × 2 5D-to-11D horizontal lifts = **1,224** output rows.
- Phase of the external forcing is recovered and checked against the original physical 10D source; no arbitrary forcing phase is injected.
- The 11D tangent generator is the numerical Jacobian of the historical material RHS at *fixed external phase*, with a midpoint short-horizon variational exponential; this is explicitly a finite-lag approximation, not a re-integrated exact tangent flow.
- The 5D finite-lag polar rotation is projected onto the 4D orthogonal complement of the physical velocity. Its exterior-square action on the 10D bivector is then compared to the y minimum-turn transport.
- Additional gauge-covariance and half-horizon diagnostics are required to distinguish transport from frame convention or discretization.

## Artifacts

- `results/model_tangent_complement_transport.csv`: 1,224 model/minimal comparisons.
- `results/authentic_source_and_jacobian_checks.csv`: forcing recovery, 11D tangent stencil, conditioning, frame quality.
- `results/model_frame_gauge_covariance.csv`: 48 covariant gauge checks.
- `results/finite_lag_horizon_sensitivity.csv`: 12 finite-horizon sensitivity probes.
- `results/selected_tangent_transport_10d_45d.npz`: reconstructible selected directions and oriented 45D polygons.
- `results/model_tangent_complement_{original,autonomous}.png`: summary figures.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`: conclusions and limitations.

Python 3.13, NumPy, SciPy, Matplotlib. The full y server sendback, including nested historical Python provenance and physical traces, is incorporated byte for byte.
