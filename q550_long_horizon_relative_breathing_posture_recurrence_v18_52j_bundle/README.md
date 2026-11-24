# Q550 — v18.52j: Long-Horizon Relative Breathing & Posture-Recurrence Audit

## Purpose

Authenticated continuation of the **same nonlinear 10D model**, now from the
v18.52i endpoint at **global xi = 3.2** to default **xi = 6.4**. No force, model
parameters, midpoint kappa, forcing phase or physical five-dimensional metric
are changed. Reference and kick integrate separately in global xi, and their
physical 5D relative vectors are compared at matched material time tau.

The archived i audit yielded 16 local D minima, 8 local maxima and near-identical
normalized waveforms across 8 kick cases. This is motivation for a longer
window, **not** a confirmed period, material memory, or closed orbit.

## Server command

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_long_horizon_relative_breathing_posture_recurrence_v18_52j_bundle.zip
cd q550_long_horizon_relative_breathing_posture_recurrence_v18_52j_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --xi-end 6.4 --samples 16385 --refine \
  2>&1 | tee ../v18_52j_server.log
```

The previously used external v18.51s server archive must remain at its
historical location, with SHA-256:
`553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`.
To locate it elsewhere, pass `--source-s-path /absolute/path/to/archive.zip`.

To validate inputs without nonlinear integration:

```bash
python3.13 scripts/run_audit.py --mode preflight
```

For an intermediate probe without modifying model physics, use
`--xi-end 4.8 --samples 8193 --refine`; after inspection, rerun with the full
horizon. The main default horizon is 6.4, with safety cap 9.6.

## Genealogy and products

- Full authenticated i server sendback under `inputs/`, including h, g and
  deeper f/e/d lineage and original scripts.
- Original byte-pinned i generator under `reference/`.
- `scripts/run_audit.py`: actual new nonlinear generator, material-time
  comparison, event search and physical 5D diagnostics.
- Full reference and kick 11D state arrays (10D+material clock) for nominal and
  refined flows, along with 5D `dQ,dV,da` arrays in NPZ.
- CSV event/ledger/summary, JSON full report, two diagnostic figures and output
  hash manifest. A complete sendback ZIP is created next to the bundle folder.

## Primary questions

1. Do additional sign-changing zeros of G=D' produce new minima/maxima of D?
2. Are the intervals between successive extrema approximately stable as the
   material-time horizon grows? Do not infer exact periodicity from coincidence.
3. Does the *fixed h-seam* angular anchor in 5D progress even when D returns?
4. At successive D minima, how do the **actual 5D vectors** differ from the
   second i-minimum in angle, norm ratio and Euclidean residual?
5. Does a near-equality of radii nevertheless retain a substantial posture
   difference, phase mismatch or hidden 10D state mismatch?

## Scientific guardrails

All zeros are numerically bracketed in a finite grid and independently refined
by DOP853; they are **not interval-certified exhaustive roots**. No 5D posture,
10D state, drive-clock/phase, orbit closure, material memory or periodicity is
certified by a repeated minimum of the scalar D. The angle to the fixed h-seam
anchor is carried unchanged across i and j; the j-seam angle is also exported.
The inherited i seam mismatch was at most ~4.7e-5 in relative D, so the audit
checks but does not pretend exact scalar seam coincidence.

## Requirements

Python 3.13, NumPy, SciPy, Matplotlib; run with OpenBLAS threads set to 1.
Only lightweight authenticity checks, synthetic tests and mock restart have
been run locally. Nonlinear server results remain **unknown** until sendback.
