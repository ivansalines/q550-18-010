# Q550 v18.51u — Scientific decision rules (specified *before* server run)

A full audit is a **PASS** only if all of the following are true for **both** authentic midpoint seed families:

1. Complete t-parent CRC, source SHA256 and output SHA256 are consistent, and nested s/r physical archives can be read; parent scientific status is `full`/PASS.
2. The tangent generator and parameter sources have the expected analytic-model provenance and matching 257/513-node grids.
3. All `u > xi` entries in the constructed two-time causal kernel are zero, and current-time instantaneous positional source feedthrough is below `1e-9` absolute.
4. Signed channel recomposition error is smaller than `5e-11` in radial/azimuthal sensitivities.
5. Reintegrating the kernel to cumulative responses matches v18.51t's **separately forced** channel IVP at each saved node within `2e-7` absolute.
6. Independently integrated tangent impulse restarts agree with Green kernel response at the final observation to `2e-7` absolute.
7. Coarsening the 257-node generator/sources to 129 interpolation nodes alters terminal response by less than `2e-5` absolute.
8. Coarsening the final source-time Simpson quadrature from 513 to 257 nodes changes terminal response by less than `2e-5` absolute.
9. Source-age weights are finite and the 1% radial compensation threshold is actually observed; if not, the audit reports a failure of the proposed marker rather than inventing a birth point.
10. The fundamental propagation remains sufficiently conditioned (`cond(Phi)<1e11`).

The separate post-run results test additionally checks that the eight terminal 0.025-width source windows exactly reconstruct each signed terminal response for all channels and both projection types.

### Classifications

- **Full PASS:** `TWO_TIME_CAUSAL_GREEN_NUMERICAL_DIAGNOSTIC_NOT_MATERIAL_MEMORY_CERTIFICATE`.
- **Preflight:** parent authentication only; no source-time kernel or birth result has been computed.
- **Failed full run:** preserve error output and do not describe the result as a certified physical birth or cross-channel material causal mechanism.

These tolerances are numerical audit thresholds, *not* validated interval error bounds or assertions of exact physical causality.
