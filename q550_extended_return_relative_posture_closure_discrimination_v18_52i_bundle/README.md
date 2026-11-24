# Q550 — v18.52i: Extended Return, Relative Posture & Closure Discrimination Audit

## Purpose

Continue the **same authentic 10-dimensional nonlinear physical model** after the
SHA-authenticated v18.52h endpoint at global `xi=1.6`. Default final horizon:
`xi=3.2` (double the previous total reach); the direct continuation spans 1.6
units of global xi. Keep the historical fixed midpoint kappa, forcing, physical
5D Cartesian embedding, and independently integrated material clocks.

The historical v18.52h observations (8/8 positive-to-negative crossings of
`G=D'`, 8/8 local maxima of `D`, and positive `C=D''` by `xi=1.6`) are
**inputs**, not predictions of what happens after the seam.

## Reproducible server command

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_extended_return_relative_posture_closure_discrimination_v18_52i_bundle.zip
cd q550_extended_return_relative_posture_closure_discrimination_v18_52i_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --xi-end 3.2 --samples 8193 --refine \
  2>&1 | tee ../v18_52i_server.log
```

A SHA-authenticated `v18.51s` external sendback is still required under
`/root/elab1/hysteresis/q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`, as for g and h.
If stored elsewhere pass `--source-s-path /absolute/path/to/that_file.zip`.

To inspect the parent without nonlinear integration:

```bash
python3.13 scripts/run_audit.py --mode preflight
```

If an extended integration fails or proves too expensive, do not modify the
physics. A shorter, clearly declared horizon can be explored with
`--xi-end 2.4 --samples 4097 --refine`. Rerun the entire command when changing
the horizon; output files are overwritten by each run.

## Source and output

- `inputs/...v18_52h_bundle_server_sendback.zip`: complete audited parent archive,
  including its g archive and deeper executable genealogy.
- `reference/parent_v18_52h_run_audit.py`: **byte-identical** real h generator,
  verified against the embedded parent ZIP.
- `scripts/run_audit.py`: exact v18.52i numerical generator.
- `tests/`: synthetic sign-crossing, physical-geometry and global-clock restart
  tests, followed by a post-run result integrity test.
- `results/extended_10d_states.npz`: actual reference and kick trajectories,
  coarse and independent refined flows, plus the global xi grid.
- `results/extended_5d_vectors_and_traces.npz`: all matched-material-time
  `dQ`, `dV`, `da` physical vectors and the radial/tangential split.
- `results/extended_ledger.csv`, `extended_events.csv`, `extended_summary.csv`, `report.json`:
  bracketed roots, decisions and limitations.
- `results/extended_return_relative_posture.png` and
  `projected_relative_position_5d.png`: diagnostic illustrations; the projection
  **does not** certify complete configuration closure.
- `SOURCE_MANIFEST.json`, `OUTPUT_MANIFEST.json`: reproducibility checksums.
- A complete sendback ZIP is written next to the bundle folder after success.

## Classification caution

A new zero of `G` from negative to positive establishes a sampled-bracketed
local **minimum of squared 5D relative position separation**, not a physical
orbit closure. Equal `||dQ||` does not mean equal 5D direction, equal relative
velocity, equal hidden 10D posture, or equal drive phase. A full material
memory/closure statement requires additional independent tests.

The code stores physically meaningful 5D vector components so that a future
audit can test multi-dimensional orientation and state recurrence explicitly.

## Environment

Python 3.13, NumPy, SciPy and Matplotlib. Check `requirements.txt`. Use
`OPENBLAS_NUM_THREADS=1` to avoid incidental multithreading overhead. Only
lightweight preflight and synthetic validation are intended on the local client.

## Patch 1 — export-only name mismatch (2026-11-07)

The original v18.52i source completed all eight physical continuations but failed
at `extended_ledger.csv` with `KeyError: 'angle_from_h_seam_deg'`. The in-memory
angle was already present under `angle_from_seam_deg`, consistent with the plot.
Patch 1 retains the CSV column `angle_from_h_seam_deg`, adds an explicit read
alias, and exercises it in the preflight regression test. No physical equation,
parameter, phase, solver, sample count, or decision rule has changed.

**Recovery:** the two large NPZ products were written before the failure.
They are inspectable, but the report, event ledger and refinement certification
were not finished. To obtain a complete authenticated sendback, rerun
`bash run_on_server.sh --xi-end 3.2 --samples 8193 --refine` with this patch.
The existing NPZ results may be backed up; this version does not silently
convert an incomplete run into a full certificate.
