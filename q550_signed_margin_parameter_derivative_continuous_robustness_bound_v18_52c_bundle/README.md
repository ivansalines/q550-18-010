# Q550 v18.52c — Signed-Margin Parameter Derivative & Continuous-Robustness Bound Audit

## Purpose

Turn the 26 authenticated v18.52b physical branch samples into a detailed **parameter-derivative audit** and a **conditional two-dimensional positivity bound**, while testing each kappa gap with a genuinely new nonlinear+tangent midpoint replay. This package deliberately does **not** claim an unconditional continuous robustness certificate.

## Inputs

- Embedded, exact v18.52b server sendback (4.6 MB), including the physical branch array, source scripts, and its 52a ancestry.
- The **original exact** v18.51s server sendback is required for the physical midpoint replay. The light bundle expects this ZIP in the parent hysteresis workspace or under `inputs/`; the full bundle embeds it. The expected SHA-256 is checked before use.

## Server command

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_signed_margin_parameter_derivative_continuous_robustness_bound_v18_52c_bundle_source_only_bundle.zip
cd q550_signed_margin_parameter_derivative_continuous_robustness_bound_v18_52c_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52c_server.log
```

Default: **24 true physical midpoint branches** (`--probes all`): 12 parameter gaps in each of two seeds. The historical nonlinear model and 10D forced tangent are recovered from the authenticated genealogy and rerun. This work belongs on the user server.

Optional `bash run_on_server.sh --probes outer` tests only the outer two parameter gaps per seed (4 physical branches). `--probes none` skips all new physical replays and is never considered a full result. For a lightweight read-only inspection, use `python3.13 scripts/run_audit.py --preflight-only`; this does not read the large v18.51s input or run any physical integration.

## Results

- `results/parameter_gap_derivative_ledger.csv`: real parameter gap widths, finite secants and surrogate midpoints.
- `results/LOCAL_PREFLIGHT.json`: the 26-branch historical analysis and explicit conditional derivative assumptions.
- `results/physical_midpoint_probe_ledger.csv`: new independent nonlinear/tangent midpoint checks.
- `results/physical_midpoint_margin_curves.npz`: signed margin as a function of xi for each new midpoint.
- `results/signed_margin_parameter_derivative_conditional_bound.png`: diagnostic figure explicitly labeled conditional.
- `results/physical_midpoint_interpolation_defects.png`: mismatch of actual physical midpoint curves vs the PCHIP surrogate.
- `results/report.json` and `results/SCIENTIFIC_REPORT_SERVER.md`.
- `SOURCE_MANIFEST.json` and fully enumerated `SHA256SUMS_OUTPUT.json`.
- The compact server sendback, created beside this bundle folder.

## Scientific honesty

The parameter **secants** and PCHIP derivatives are not automatically rigorous upper bounds for derivatives of the nonlinear model. Positive midpoint tests do not eliminate unsampled sign changes. The putative derivative inflation factor defaults to 64 and is clearly marked **ASSUMED**, not certified.

Historical v18.52b's output manifest omitted five generated files; this package hashes all five directly and records their exact SHA-256 values. The v18.52c output manifest enumerates all its generated results.

All source code, comments and scientific documentation are in English. The ZIP and its internal root share the full descriptive audit name.
