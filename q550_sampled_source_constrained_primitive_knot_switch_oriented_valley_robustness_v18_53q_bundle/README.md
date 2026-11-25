# Q550 — v18.53q: Sampled-Source-Constrained Primitive, Knot-Switch & Oriented Valley Robustness Audit

This audit continues v18.53p and tests the numerical representation of the full-10D residual valley. It uses the authenticated v18.53p server sendback, including its exact v18.53o physical ancestry. **No new Q550 IVPs are performed**. Its goal is to distinguish a nonsmooth primitive/argmin branch from a persistent 10D orientation change.

Scientific inputs: 2,241 archived clock samples for each of six physical cases, with cumulative 10D primitive C(xi), sampled bivector source B_10D, and clock-rate factor. All source files are SHA-256 and CRC authenticated. The script extracts the nested original-o numerical arrays from the pinned parent-p ZIP, without rewriting the originals.

## Checks

1. Authenticate p and nested o, including complete source/output manifest hashes. Confirm no source was regenerated.
2. Compare three C reconstructions from precisely the same samples: componentwise PCHIP(C), cubic Hermite(C,C'=sampled physical source), and integral of piecewise-linear sampled source (source-integrated quadratic, anchored at C(xi0)). The latter two are **interpolants**, not newly integrated physical trajectories. A source-integrated curve need not match all stored C knots because the cumulative export can have quadrature differences; report that drift rather than forcing agreement.
3. At each origin, search for **all bracketed stationary widths** by the analytic derivative of ||N||^2; select global bounded minimum among candidates and report ties, active bounds, source curvature, and near-clock knots.
4. Compare the 10D unit-direction path and its 45D oriented wedge integral. Use both direct geodesic polygonal length and chain-rule tangent length. For interior minima, derive w'(origin) by implicit differentiation of N·g_end=0; flag the ill-conditioned cases where the denominator nearly vanishes.
5. Refine 141→281→561 origins; compare all six physical case histories. No result in xi>25.6 is called independently measured material return.

## Run on server

```bash
cd /root/elab1/hysteresis
B=q550_sampled_source_constrained_primitive_knot_switch_oriented_valley_robustness_v18_53q_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 bash "$B/scripts/run_on_server.sh" 2>&1 | tee v18_53q_server.log
```

Optional more width stationary-point brackets: add `--deep-search`. The run generates an archive with full descriptive name ending `_server_sendback.zip`. To test locally without the full postprocessing, execute `python3 scripts/run_audit.py --preflight` and `python3 tests/test_synthetic.py`.

## Readout

`results/report.json` and `results/SCIENTIFIC_REPORT_SERVER.md` classify the numerical outcome; `results/convergence_summary.csv`, `results/valley_samples.csv`, `results/cross_interpolant_comparison.csv`, `results/primitive_source_consistency.csv`, and `results/oriented_forms.npz` retain the measurements. The parent ZIP is unchanged and preserved inside the distribution/sendback. All scientific texts, comments and reports are in English.

## Limits

The two source-constrained reconstructions are constrained by **sampled** physical-source values only, not dense physical replay between knots. The Euclidean 10D bivector norm is not the full 11D state metric. Neither near-return, curvature, nor 45D oriented transport establishes phase locking, material holonomy, or periodic closure. A truly independent test of source-density effects requires new server-side physical integration at a finer clock grid, after this numerical discrimination.
