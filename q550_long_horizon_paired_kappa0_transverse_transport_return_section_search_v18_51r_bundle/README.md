# Q550 — v18.51r

**Long-Horizon Paired-κ₀ Transverse Transport & Return-Section Search Audit**

## Purpose

Continue the physical nonlinear Q550 two-κ₀ comparison beyond the v18.51q local horizon without conflating the two independently prepared seeds, without treating a local phase projection as orbit shadowing, and without presupposing a material return. The horizon is fixed in advance at `0 <= xi <= 0.2`, four times the parent audit's span. A failed/incomplete physical IVP is a negative audit outcome and must be reported, not replaced by extrapolation.

## Historical source and experimental design

The nested parent is the complete authenticated `q550_paired_kappa0_shadowing_phase_lag_error_budget_v18_51q_bundle_server_sendback.zip` with executable genealogy to the original 10D driven physical solver. It is verified by both source and output SHA-256 manifests. The authentic, original, autonomous, and midpoint κ₀ values are read from the parent model, never fitted.

Both independent initial seed preparations are integrated at the same three κ₀ choices. Each pair is replayed at three DOP853 numerical resolutions (`standard`, `reference`, `tight`), 18 distinct integrations in total. These three solver settings form **numerical discrepancy probes**, not interval bounds. Each branch saves 513 event-relative sections, full Cartesian position/velocity and 10D material states, clocks and polar phases.

## Operational observables

- Same-seed cross-parameter Cartesian `delta Q`, `delta V` and all-tangent state differences.
- The local along-orbit time-equivalent lag and the orthogonal residual as separate quantities.
- A *single* constant time lag fitted on the predeclared `[0.05,0.195]` window, independently from local per-section shifts.
- A search for zeros of `dx/dxi` and for recurrence of the initial projected `x` coordinate, with full-state distances for any detected equal-projection event. A projected return is not full-state closure.
- Difference between three solver settings and the overlap at `xi=0.05` with the authenticated v18.51q physical-state endpoint.

See `MATHEMATICAL_PROTOCOL.md` and `SCIENTIFIC_DECISION_RULES.md`. Scientific documents and code comments are in English.

## Run

```bash
cd /root/elab1/hysteresis/q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_51r_server.log
```

The runner validates the nested parent, performs independent synthetic tests, replays the 10D model, independently checks result identities, writes SHA-256 manifests and packages a same-root `_server_sendback.zip`. If extended IVPs fail to complete, retain and share the full log; do **not** claim a successful r audit.

## Caution

`R` in this version label means the **next audit letter**, not a validated radial coordinate or a new material surface. The historical toroidal coordinate `R`, the physical scalar Q=550, and the Cartesian immersion `Q=(Q1,...,Q5)` are distinct objects.
