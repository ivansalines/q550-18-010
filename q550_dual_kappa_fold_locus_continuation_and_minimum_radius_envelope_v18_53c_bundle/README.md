# Q550 — v18.53c
## Dual-κ₀ Fold-Locus Continuation & Minimum-Radius Envelope Audit

**Purpose.** Continue the genuine v18.53b `Δ = π` relative-equilibrium saddle-node across the *two distinct* historical κ₀ prescriptions and ask whether the critical **smaller equilibrium radius** has an interior minimum/maximum between them. The common fixed charge is the authenticated v18.53a/v18.53b value, **J = −1.441016445824766**. The signed drive `s` is held constant *within each equilibrium equation* but its fold location `s_fold(κ₀)` changes during continuation.

**Historical values:**

- κ₀ original `1/6 = 0.16666666666666666`.
- κ₀ autonomous `0.1665996104` (historical rounded value, not a new fitted coefficient).
- κ₀ common v18.52l integration `0.1666331385421536`; v18.53b authenticated the fold there, at `s_fold≈0.9268428336365817`.

The original and autonomous *seed offsets* at the common κ₀ were gauge/phase equivalent: they were not two independent physical solutions at the common anchor. Here original and autonomous refer to **different κ₀ parameter prescriptions**, and their fold equilibria generally occur at **different signed drives**.

### Run on the user's server

From `/root/elab1/hysteresis`:

```bash
unzip q550_dual_kappa_fold_locus_continuation_and_minimum_radius_envelope_v18_53c_bundle.zip
cd q550_dual_kappa_fold_locus_continuation_and_minimum_radius_envelope_v18_53c_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --require-parent 2>&1 | tee ../v18_53c_server.log
```

If the SHA-pinned v18.51s ancestor is stored elsewhere:

```bash
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  --require-parent \
  --parent-v1851s /absolute/path/to/q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip
```

The original v18.51s ZIP is external and SHA-pinned, not duplicated into this small bundle. The v18.53b authenticated scientific report, its refined fold and the v18.53a report are embedded with SHA-256 checks. The actual **v18.53b algorithm** and earlier source files are bundled, not just cited.

A stand-alone diagnostic, without source authentication, is available via `bash run_on_server.sh --offline-only`. It must report `OFFLINE_FOLD_LOCUS_SCREENING__ORIGINAL_ENGINE_NOT_VERIFIED` and should **never** be described as authentic Q550 field certification.

### Numerical workflow

1. Verify all bundled source SHA-256 hashes, validate the analytic fold Jacobian against finite differences and run independent archived RHS probes.
2. Recheck the authentic v18.53b fold and solve the 4D augmented algebraic system `F(a;J,κ,s)=0, det ∂_aF=0`, with `a=(log r1,log r2,ω1)` and `J` fixed.
3. Continue **102 high-resolution sampled points** between the two κ₀ prescriptions (101 grid points plus central authenticated κ₀) using the implicit derivative as predictor.
4. Optionally continue **64 wider sampled points** over `0.085≤κ₀≤0.25` to see the larger fold geometry. Do not misclassify a corrector failure as a fold termination.
5. At every accepted point, check the nondegenerate radial fold, one real near-zero reduced mode and stability of the other seven modes. With authenticated parent, check **original 11D RHS and reduced 8D Jacobian** at every root.
6. Report `min(r1_fold,r2_fold)` at both endpoints and across the sampled interior; compute `dr_min/dκ₀` using implicit differentiation. Do not claim a rigorous continuum envelope from finitely many samples.
7. Solve the original and autonomous κ₀ equilibria **also at common fixed drive `s=-1`**, a different comparison from the two critical fold radii.

### Files returned

- `results/report.json` — verdict, historic provenance, endpoint radii and signed drives, observed envelope, derivative tests, same-drive comparison and explicitly limited claims.
- `results/fold_locus_continuation.csv` — high-resolution and extended sampled curves, radii, signed-drive fold, implicit slopes, native checks and critical spectral data.
- `results/fold_locus_selected_states.json` — complete critical core states and local fold diagnostics.
- `results/continuation_stops_not_proven_fold_termination.json` — any corrector/continuation stops.
- `results/dual_kappa_fold_radius_envelope.png`, `results/dual_kappa_minimum_radius_zoom.png` — visual comparisons.
- `results/SCIENTIFIC_REPORT_SERVER.md` — interpretation and exclusions.
- The full source, `SOURCE_MANIFEST.json`, `OUTPUT_MANIFEST.json`, tests, and original SHA-pinned historical interfaces.

**Not tested:** historical forward–return trajectory, two different dynamical minimum radii, global radial noncollapse, global fold-locus topology, and interval-certified monotonicity. All source files and mathematical comments are English by design.
