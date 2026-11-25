# Q550 — v18.53w
## Fixed-Center Metric–Width Separation & Moving-Minimum Oriented Geometry Audit

### Purpose
Disentangle the previously reported ~6.56 degree 45D oriented-polygon change from a genuinely moving window center. **In v18.53v, the oriented polygon was always sampled on an origin grid centered on the canonical kinetic minimum (`oc`)**; only the metric and the optimized width at each origin changed. The metric-dependent *global minimum* moved, but that displacement was **not part of the v18.53v 45D polygon**.

This bundle introduces an explicit **three-factor factorial audit**:
1. `center`: origin grid centered on canonical min or metric-conditioned min.
2. `width`: width minimization with the canonical or the metric-conditioned 10D quadratic norm.
3. `metric`: square-root Euclidean identification under the canonical or metric-conditioned induced exterior-square metric.

For each of the 2^3 combinations, the same source-dependent 10D accumulated bivector is evaluated on the corresponding window. The four curves corresponding to choices of center and width are compared using both metrics. The historical v18.53v `center=0,width=1,metric=1` polygon is replayed **from its archived optimized width CSV**, not newly fitted.

### Scientific meaning
The **45-component oriented polygon** is calculated from unit 10D vectors following the archieved `oriented_polygon()` prescription. The 45D polygon vectors are all identified through the positive definite symmetric matrix square root, one **convention-dependent** common Euclidean comparison frame.

The three factor contributions are permutation-averaged Shapley *vector differences*: they sum exactly to the change in the oriented polygon. Their **norms and associated angles do not add**. Signed projections onto the total displacement can be negative or greater than one; they are not physical energy percentages. This experiment does not identify a material angular momentum, causal effect, or physical holonomy.

A second **2-factor** audit applies the same distinction to the **10D residual at the two metric-conditioned minima**, using `sqrt(H) N(p)` (center shift vs metric). It also produces a scalar Shapley decomposition of the *change in the norm*, with all cross-evaluations performed.

### Data and reproducibility
The original, unaltered v18.53v server sendback is in `inputs/v18_53v_authenticated_server_sendback.zip`. It embeds the authenticated v18.53u and v18.53t sendbacks. The original u 5D kinetic metric helper is copied byte-identically and authenticated against the ancestor's SHA-256. All source files, input ZIP, mathematical protocols, checks, outputs and manifests are bundled.

**No physical IVP is executed**; the 8961-point authentic 21D traces are reused. The full run performs 204 metric continuation nodes, with 141 curve samples per node. This is more expensive post-processing than the smoke tests but is intended for the user's server.

### Execution (Python 3.13)
```bash
cd /root/elab1/hysteresis
B=q550_fixed_center_metric_width_separation_moving_minimum_oriented_transport_v18_53w_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
  bash "$B/scripts/run_on_server.sh" --require-parent \
  2>&1 | tee v18_53w_server.log
```

Other modes: `bash "$B/scripts/run_on_server.sh" --preflight-only` or `--smoke-only`.

Outputs: `metric_path_factorial_45d.csv`, `metric_path_factorial_vectors.npz`, `metric_path_family_summary.csv`, two figures, `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`, and `...v18_53w_bundle_server_sendback.zip` next to the extracted bundle.

### Limitations
- The kinetic canonical metric is identity in the historic mass-whitened Cartesian embedding; the velocity-occupation metrics are *noncanonical diagnostic choices*.
- The adapted center and its optimum width derive from the same interpolated trajectory, not new trajectories.
- Far endpoints of the three-window path lie beyond the xi=25.6 preformal clock boundary: the analysis is in **formal clock continuation**.
- This audit can separate the observer/representation contribution and path-selection contribution in the chosen convention. It does **not** establish material closure, memory, or holonomy.
