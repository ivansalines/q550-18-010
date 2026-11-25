# Q550 — v18.53j: Three-Window Residual-Bivector Near-Closure, Window-Shift Robustness & Temporal Holonomy Audit

Portable, SHA-anchored, executable audit for the user's Linux Q550 server. **All code and documentation are in English.** The large original ancestor and native historical files are authenticated externally; they are not stuffed into this small ZIP.

## Motivation

The authenticated v18.53i server run completed **26/26** independent Q550 IVPs to ξ=51.2 and observed three 5D *time-integrated residual torque bivectors* with approximately equal magnitudes and mutual bivector-space angles near **120°**, yielding a small additive net. The original seed / shared κ₀ baseline had ≈93.7% compensation **among the three time windows**, different from the ≈98.4% whole-interval cancellation **between the radial and phase torque integrals**. These two metrics must never be conflated.

The present experiment asks whether the near-triplet geometry is robust to predefined *whole-window translations*, *window-duration changes*, *internal boundary jitters*, *historical vs autonomous κ₀*, and initial historical seed. It explicitly tests whether the alleged temporal holonomy can be defined independent of arbitrary Schur-frame gauges. We make **no claim** that an integrated torque nearly returning to zero is a physical orbit returning to its state.

## Predeclared experimental design

- Starting states: both original and autonomous native v18.52l historical seeds at **ξ=12.8**. Each seed supports **three** genuine IVPs: shared κ₀ = 0.1666331385421536, original κ₀ = 1/6, and autonomous κ₀ ≈ 0.16659961041764054. The *same external phase forcing* applies to the three κ variants within a seed.
- The six primary trajectories run to **ξ=56.0**; **four higher-accuracy original/autonomous-κ₀ refinements** serve as integrator checks. Expected total: **10 completed IVPs**, not 26.
- Original native archive: **ξ=12.8–25.6**, first replay checked against native SHA. The v18.53i baseline **ξ=12.8–51.2** is independently reconciled to exact 53i outputs per segment. The range **ξ>25.6**, including ξ=51.2–56.0, is an extrapolation of the authenticated generator/clock protocol, NOT an independently observed historical continuation.
- First scan: whole interval origin `ξ₀∈{12.8,13.6,14.4,15.2}` and equal window width `w∈{9.6,10.8,12.0,12.8,13.6}`. Three intervals: `[ξ₀,ξ₀+w]`, `[ξ₀+w,ξ₀+2w]`, `[ξ₀+2w,ξ₀+3w]`; 20 triplets × six dynamics = **120** independent windows-scoring comparisons.
- Second scan: fixed entire endpoints `[12.8,51.2]`, but first cut at `25.6+δ₁` and second at `38.4+δ₂`, with `δ₁,δ₂∈{-2,-1,0,1,2}`. 25 triplets × six dynamics = **150** diagnostic comparisons. Their NET residual is *necessarily identical* by algebraic telescoping; only the individual intensities, orientations, angles and normalized score can reveal cut sensitivity. Do NOT count fixed-endpoint net invariance as physical robustness.
- In all **270** triplets, compute the complete 5×5 bivectors, their Frobenius norms, pairwise angles, three-window near-closure `1−||R₁+R₂+R₃||/(||R₁||+||R₂||+||R₃||)`, magnitude CV, radial-vs-phase cancellation, and numerical diagnostics.
- At each triplet, optionally construct a *closed geodesic triangle of principal 2-plane subspaces* in the Grassmannian, transporting Schur frames by a minimal-overlap polar alignment. Return angle is an **observer-space GEODESIC POLYGON diagnostic**, not a naturally closed physical time-evolution and not a Noether charge. Report undefined when principal strengths become degenerate or adjacent principal subspaces are nearly perpendicular. Do not silently supply the missing physical closure.

## Acceptance / failure criteria

1. Source references v18.53g, h fix1, and i all pass pinned SHA256 and previous output-manifest authentication; v18.53i report verifies 26 genuine IVPs. The original v18.51s and v18.52l SHA256 are verified before entering the scientific branch.
2. Reconstructed authentic material-time clock `dτ/dξ=ω/sigmoid(z)` passes archived rate validation.
3. Every of the six new trajectories has source agreement below 2e-10, conserved counter-rotation charge drift below 5e-9, positive material time, and real finite states.
4. Shared-κ₀ trajectory reproduces the original native prefix at scaled max error below 2e-4 and four refined integrations agree with standard integrations below 3e-7.
5. Every *historical v18.53i three-window reference segment* is reconstituted within a 5e-4 absolute residual norm tolerance and 2e-4 absolute cancellation tolerance; these are conservative, declared before running. Independent off-grid material-time quadrature is also checked.
6. The numerical near-closure threshold is **90%**, equiangular tolerance ±8° around 120°, and equal-strength CV threshold 10%. Results *below* thresholds are not rejected as failures of data authenticity; they are **scientific falsifications** of an extended persistence interpretation.
7. If large ancestors are absent, the audit writes a clear **partial** report and sendback. With `--require-parent --require-native --require-dynamic`, it exits nonzero. If the code fails after some IVPs, the report keeps the actual count and audit status.

## Server execution

```bash
cd /root/elab1/hysteresis
B=q550_three_window_residual_bivector_near_closure_window_shift_robustness_temporal_holonomy_v18_53j_bundle
unzip "${B}.zip"
cd "$B"
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
 bash run_on_server.sh --require-parent --require-native --require-dynamic \
 2>&1 | tee ../v18_53j_server.log
```

Explicit ancestor flags are also supported: `--parent-v1851s /path/to/authentic/parent.zip` and `--native-npz /path/to/extended_horizon_10d_states.npz`. The trusted loader searches the parent directory recursively when omitted.

## Outputs

- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`: scientific classification, actual IVP count, authentic source checks, objective thresholds, and explicit limitations.
- `results/window_shift_width_triplets.csv`: 120 true whole-interval changes, with angles, strengths, near-closure and polygon transport.
- `results/fixed_endpoint_internal_boundary_jitter.csv`: 150 boundary-jitter comparisons at identical overall endpoints; explicitly *no physical inference from constant net*.
- `results/all_residual_bivector_components.csv`: all 810 individual 5D skew residual matrices in upper-triangle coordinates plus separate source magnitudes. Enough to independently reconstruct the near-closure and principal-plane geometry.
- `results/v1853i_three_window_reconciliation.csv`: all 18 per-segment baseline comparisons versus SHA-pinned v18.53i.
- `results/integrated_ivp_replay.csv`, `integrator_refinements.json`, `v1853g_clock_reconciliations.json`, `offgrid_quadrature_reconciliation.csv`: numerical and genealogy tests.
- `results/baseline_three_window_geometry.json`: each seed / κ baseline measurement.
- Four optional `.png` plots show whole-interval shift/width near-closure and fixed-endpoint angular distortion for the two seed families.
- `SOURCE_MANIFEST.json`, `OUTPUT_MANIFEST.json`, and parent-level `..._server_sendback.zip`: reproducibility.

The complete v18.53i source `run_audit.py` and all prior h/g/f physical generator support and source references remain **inside the ZIP**. No external web access, background process, or undocumented modeling change is needed. In `MATHEMATICAL_PROTOCOL.md` we explicitly distinguish the (Euclidean Frobenius) bivector-space 120° from a three-dimensional spatial angle.
