# Q550 — v18.53m-fix1
## Two-Parameter Fixed-Frame Near-Null Locus, Conditioned Jacobian & Boundary-Source Attribution Audit

Fix1 repairs the numerical optimizer search domain only. The server archive v18.53m failed after its first Q550 IVP with `ENDPOINT_OUTSIDE_TRACE` because origin bounds began at xi=12.6, but the native trace begins at xi=12.8. That authenticated failed return is preserved in `provenance/`; this package keeps all scientific tests and original physical inputs unchanged. See `FIX1_CHANGELOG.md`.

This audit follows authenticated Q550 v18.53l. The question is whether **both independent coordinates of the integrated net bivector projected on one fixed rank-two 10D bivector-space plane** can vanish simultaneously as the *origin* `o` and the *common three-window width* `w` vary. We do not mistake a small norm on a one-dimensional width scan for a 2D root; we also refuse to confuse any projected zero with a full 10D physical closure.

### Run on research server

```bash
cd /root/elab1/hysteresis
B=q550_two_parameter_fixed_frame_near_null_locus_conditioned_jacobian_boundary_source_attribution_v18_53m_fix1_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
bash "$B/scripts/run_on_server.sh" \
  --require-parent --require-native --require-dynamic \
  2>&1 | tee v18_53m_server.log
```

The script discovers the original historical v18.51s parent and authenticated native 10D dataset under the research workspace, as in v18.53l. Or pass their exact locations explicitly:

```bash
PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" \
  --parent-v1851s /ABS/PATH/authentic_parent.zip \
  --native-npz /ABS/PATH/authentic_native_10D.npz \
  --require-parent --require-native --require-dynamic
```

**Do not interpret preflight as physics.** To verify archive, mathematical unit tests and ancestor integrity without integrating:

```bash
PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" --preflight
```

It produces an explicitly labeled `PREFLIGHT_PASS__NO_NEW_Q550_INTEGRATIONS` report. The main physical run requires both historical large files and checks their exact SHA-256 pins before doing *any* new Q550 integration. No synthetic data or CSV reuse substitutes for them.

### Research geometry

A skew 5×5 residual bivector is embedded isometrically in ten coordinates by scaling upper-triangular entries by sqrt(2). The original source torque contributions, original material clock, and three-window trapezoidal primitive are taken from the actually integrated Q550 model.

Choose a **single 10×2 orthonormal basis** `E` from the SVD of the three original integrated bivectors at `(o_anchor,w_anchor)`, keep it fixed for *both* variable parameters, and solve

`F_E(o,w) = E.T [C(o+3w)-C(o)] = (0,0)`.

This basis depends on the chosen observer anchor; to detect that dependence, repeat with *three independently fixed bases* anchored at `o=13.6,14.4,15.2`, each at `w=12.8`. No SVD basis is reoriented to follow the root search. The `SO(2)` internal re-basing test confirms that changing the axes **within the same fixed plane** does not change norms, root positions, or singular values.

Main grid: 25 origins `o=12.8..15.2`, 29 widths `w=12.74..12.88`, all evaluated for two historical seeds × three κ₀ rules × three fixed anchor frames. Search bounds are restricted to actual supported trace: `o=12.8..15.4`, `w=12.70..12.91`. Each group uses multi-start bounded least-squares followed by uniqueness consolidation and complete conditioning measures. Width and origin are **postprocessing of six integrated traces**, not separate dynamic runs.

### Genuine server execution and validation

- 6 original Q550 11D primary integrations (two seeds × three κ₀ rules), plus 4 tightened repeats.
- Reuse the source-staging methods and native reconciliation from the **full pinned v18.53l sendback** embedded byte-identically in `inputs/`; it includes the k/j and older numerical-code genealogy.
- Verify both source and output manifests of l, the historical l decision, its 756 real measurement rows, embedded k and j authenticity, and 54 initial ten-component regression comparisons with l.
- Repeat after every-other-xi-node re-quadrature. Record any root shift, not merely accuracy of the optimizer on an unchanged integral.
- For four tightened integrations, independently recompute the torque/clock primitive, evaluate every candidate, and re-optimize conditioned roots to measure root displacement.
- Export exact boundary-slope Jacobian, source-reconstructed Jacobian discrepancy, determinant identity, SVD singular values, condition numbers in both original and domain-normalized coordinates, distance of endpoints to the nearest interpolation knot, and complete 10D complementary residual.

Expected verified output ZIP:

`q550_two_parameter_fixed_frame_near_null_locus_conditioned_jacobian_boundary_source_attribution_v18_53m_fix1_bundle_server_sendback.zip`

### Output files

- `two_parameter_fixed_frame_grid.csv`: all 13,050 two-parameter grid samples, projected two-vector coordinates and norm.
- `multistart_optimizer_candidates.csv`: every bounded solver attempt with Jacobian conditioning.
- `deduplicated_projected_root_candidates.csv`: unique solutions, alternative-anchor comparisons, and full 10D complement.
- `boundary_source_jacobian_attribution.csv`: direct two-endpoint source identity and material-sampling checks.
- `parent_l_net_reconciliation.csv`: authenticated replay comparison to previous version.
- `decimated_quadrature_root_sensitivity.csv`: every-other-node integration perturbation.
- `tight_integrator_root_sensitivity.csv`: tight-solver perturbation on four tested dynamics.
- `anchor_basis_conditioning.csv`, `genuine_integrated_ivps.csv`, `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, 2D maps and manifests.

### Scientific scope

This experiment can identify a conditioned zero of the **selected observer's two-dimensional projection**, not a full 10D zero or material return. In particular, the net three-window integral telescopes and has no dependence on the two internal partition boundaries. A two-parameter zero of `E.T[C(e)-C(s)]` is a self-intersection of the 2D projected cumulative-source curve `E.T C(xi)` at endpoints `s` and `e`; it may be a choice-of-projection effect. The data after xi=25.6 use the historical model continuation, not independently observed history.

All mathematical claims and file comments are in English. The archive is self-contained for genealogical code reproduction **provided** the two large externally pinned historical physical inputs are present on the server.
