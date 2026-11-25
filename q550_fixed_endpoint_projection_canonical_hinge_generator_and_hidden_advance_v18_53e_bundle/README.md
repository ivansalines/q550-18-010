# Q550 — v18.53e: Fixed-Endpoint Projection, Canonical Hinge Generator & Hidden-Advance Audit

**Scientific task:** test whether an observer tilted relative to the same time pair can make a relative Q550 separation appear to return while a hidden component advances. Independently reconstruct the absolute native 5D flat-kinetic hinge moment candidate and discriminate it from the exact counterrotation charge.

## Important scope distinctions

1. **Observer test:** v18.52l `dQ,dV` are perturbation-minus-reference separations. They do **not** represent one absolute orbit. The authentic v18.53d observer candidates are frozen as time pairs **before** comparing the five inclination angles, so observer comparisons have identical endpoints. An additional time pair is chosen by minimum **full 5D** distance, independent of observer. The optional analytically optimal viewing angle is explicitly post-selected and does not independently establish a physical hinge orientation.
2. **Absolute moment test:** v18.52l native **10D absolute** states are mapped with the historical v18.49h kinetic embedding into `(Q,V) in T(R^5)`. The candidate angular moment `J_hinge^(T)=(Q.n)(V.h)-(Q.h)(V.n)` uses the flat-kinetic Legendre `P=V`. It is a canonical moment of that declared **kinetic sector** but **not** a proven Noether charge or conservation law of the full driven Q550 dynamical action. The PCA-derived `h` is an observer direction, not an established physical symmetry axis.
3. **Distinct invariant:** `J_rot = r1^2 omega1 - r2^2 omega2` is the separately authenticated conserved first integral of the original Q550 angular equations. No equality with `J_hinge^(T)` is imposed.
4. **Different kappa experiments:** v18.53c authenticated fold points compare original and autonomous κ₀, each at its own fold; the v18.52l long 5D and native 10D traces use their historically shared κ₀. Never interpret the dual-κ fold variation as dynamical transport along those trajectories.

## Required external inputs (large, deliberately not copied into ZIP)

In `/root/elab1/hysteresis/q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/`:

- `extended_horizon_5d_vectors_and_traces.npz`, SHA-256 `d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`;
- `extended_horizon_10d_states.npz`, SHA-256 `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12`.

The verified small historical source copies (v18.53d output, v18.52n observer, v18.53c fold source, v18.49h kinetic embedding, and v18.52p prior hinge audit) **are included** in `reference/`, with a source manifest and previous output provenance checks.

## Run on the server

```bash
cd /root/elab1/hysteresis
unzip q550_fixed_endpoint_projection_canonical_hinge_generator_and_hidden_advance_v18_53e_bundle.zip
cd q550_fixed_endpoint_projection_canonical_hinge_generator_and_hidden_advance_v18_53e_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --require-trace --require-native \
  2>&1 | tee ../v18_53e_server.log
```

The `--base` argument can set the parent source directory. For alternate locations pass both `--trace-npz` and `--native-npz` explicitly. A missing required input produces a **partial**, non-certified report plus an uploadable sendback archive, and exits nonzero. No Q550 ODE is integrated.

## Files generated

- `results/fixed_endpoint_observer_projection.csv` — 8 relative series × 3 **fixed** time pairs × 5 predeclared inclinations plus one analytically optimized, **post-selected** viewing angle. Full 5D distance, screen distance, signed viewing depth, remaining invisible components, and velocity mismatch are recorded together.
- `results/fixed_endpoint_pair_summary.json` — time pairs and invariance against observer angle, separated by pair origin.
- `results/native_absolute_kinetic_hinge_moment.csv` — diagnostic samples from 20 verified native absolute 10D series: kinetic-plane moment, rotor Noether first integral, and finite-difference derivative diagnostics.
- `results/native_absolute_kinetic_hinge_summary.json` — ranges and variation of both moments, with full driven canonical-action caveat.
- `results/fixed_endpoint_shadow_hidden_observer.png`, `results/absolute_kinetic_hinge_vs_conserved_rotor_charge.png` — figures with explicitly different measurement scopes.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`, matching-root `..._server_sendback.zip`.

This audit is reproducible with the genuine source data and does not claim closure or an intrinsic physical hinge axis on the strength of an optimally tilted observer.
