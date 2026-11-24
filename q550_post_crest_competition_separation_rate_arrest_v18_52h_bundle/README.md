# Q550 — v18.52h: Post-Crest Competition & Separation-Rate Arrest Audit

**Purpose.** Continue the eight authentic nonlinear Q550 kick/reference trajectories **from the verified g endpoint** at global event clock `xi=0.8`, rather than re-integrating or extrapolating the g curves. Track whether the post-crest inward/kinematic competition actually arrests the squared Cartesian separation rate.

This bundle includes the **exact g server sendback** and a readable copy of the g generator. The g sendback itself contains the full f→e→d historical source and result genealogy; all files in the parent's `SOURCE_MANIFEST.json` and `OUTPUT_MANIFEST.json` are SHA-256 verified before physical calculation. Starting physical state, 10D nonlinear field, kinetic embedding, Q550 parameters, midpoint kappa and material-clock transport are unchanged. Both original and autonomous seeds, both kick signs, and both kick magnitudes are retained, with separately continued coarse and refined solutions.

## External authentic nonlinear ancestor (same as g)

The full v18.51s original engine is **not** embedded in the g archive. The same historical external dependency is required:

`q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`

Expected SHA-256: `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`.

Place it in `/root/elab1/hysteresis/` or use `--source-s-path /absolute/path/to/file.zip`. No source substitution is permitted; the run stops if it is missing or has a different hash.

## Run on the server

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_post_crest_competition_separation_rate_arrest_v18_52h_bundle.zip
cd q550_post_crest_competition_separation_rate_arrest_v18_52h_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
bash run_on_server.sh --xi-end 1.6 --samples 4097 --refine \
  2>&1 | tee ../v18_52h_server.log
```

The global `xi` horizon now extends **0.8 → 1.6**; the inherited g horizon was **0.0 → 0.8**. If the genuine nonlinear RHS or physical material clock becomes invalid, do not change the field or accept an extrapolation: retain the log and rerun a shorter horizon, e.g. `--xi-end 1.2`. A longer exploratory horizon such as `--xi-end 2.0` may be used **only after** interpreting the first result. The script refuses `xi_end > 3.0`.

The original event-clock time is preserved at restart: the drive uses `theta(event + xi)`, with `xi` initialized at `0.8` in the ODE solver. The material clock starts from the **archived**, independently integrated 11th state (not zero). A matched physical-time comparison is restricted to the **overlap** of the two branch clocks, beginning at `max(tau_base(0.8), tau_kick(0.8))` and ending at `min(tau_base(xi_end), tau_kick(xi_end))`. The few-microsecond gap at the boundary is compared to the g physical trace for continuity.

## Diagnostics and decisions

Let `q=Qkick−Qbase`, `v=Vkick−Vbase`, `a=akick−abase` in the exact Cartesian 5D kinetic representation; let `D=||q||²/2`, `G=q·v = D′`, `P=||v||²`, `I=−q·a`, `C=P−I=D″`, and `R=I/P` where P is well conditioned.

- The first g **ratio crest** is retained as a sampled, parent-verified marker (~tau 1.239), not reclassified as a new h result.
- A post-crest `C=0` crossing with **negative→positive** curvature means `R` crossed **from above to below 1**; the separation rate `G` has a local minimum, not automatically zero.
- A `G=0` crossing with **positive→negative** orientation is the *only* criterion for a genuine local maximum in `D`. A touch or unresolved event is not sufficient.
- `R` local peaks are sample-bracketed and refined against the two genuinely continued IVPs; they are not a physical-force maximum by themselves.
- A finite sampling grid cannot certify the absence of tangencies or arbitrarily fast sign changes. Independent refined IVPs check sign-changing roots, but are not rigorous interval arithmetic.

## Outputs

- `results/report.json`: parent/source provenance, per-trajectory events, checks, and conservative classification.
- `results/post_crest_ledger.csv`, `post_crest_events.csv`, `post_crest_summary.csv`: numerics and chronology.
- `results/post_crest_traces.npz`: all per-kick matched-material-time scalar traces.
- `results/post_crest_continued_10d_states.npz`: **real** 11-component (10D state + material clock) states across xi=0.8..xi_end, coarse and refined when enabled.
- `results/post_crest_competition_rate_arrest.png`: figure with the inherited ratio crest and g→h boundary indicated.
- `results/SCIENTIFIC_REPORT_SERVER.md`: scientific report, including negative findings when appropriate.
- `OUTPUT_MANIFEST.json`: output SHA-256 manifest, and the full `q550_post_crest_competition_separation_rate_arrest_v18_52h_bundle_server_sendback.zip` containing the executable generator, source, and results.

Local-only preflight (without physical continuation):

```bash
python3.13 scripts/verify_source_manifest.py
python3.13 tests/test_synthetic.py
python3.13 tests/test_mock_restart.py
python3.13 scripts/run_audit.py --mode preflight
```

**No long nonlinear integrations should be run locally; use the server.** All code comments and scientific prose in the bundle are in English.
