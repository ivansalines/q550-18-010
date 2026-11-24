# Q550 — v18.51w

**Channel-Vector Angle, Oriented-Area & Observer-Critical Compensation Audit**

Executable, genealogically traceable numerical audit of the source-to-response vectors in the second Cartesian plane, inherited from the authenticated v18.51v server sendback.

## Run on the Q550 server

```bash
cd /root/elab1/hysteresis
unzip q550_channel_vector_angle_oriented_area_observer_critical_compensation_v18_51w_bundle.zip
cd q550_channel_vector_angle_oriented_area_observer_critical_compensation_v18_51w_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_51w_server.log
```

No expensive nonlinear integrations are performed: the script authenticates the v18.51v sendback and computes a lightweight geometric postprocessing of the archived 513 sections for each seed. The parent ZIP contains the nested v18.51u/t/s/r physical-model genealogy.

## Lightweight source-only ZIP

Extract the source-only ZIP, then copy the previously obtained
`q550_physical_clock_green_transport_moving_frame_invariance_v18_51v_bundle_server_sendback.zip`
into this bundle's `inputs/` directory. The same runner then performs all verification and produces a complete sendback. The `SOURCE_MANIFEST.json` includes the full input ZIP digest and will fail if the wrong parent is supplied.

## Main outputs

- `results/channel_vector_section_ledger.csv`: 1,026 rows, including both physical clock and drive clock, signed dot/area, oriented angle and observer-critical axes.
- `results/channel_vector_geometric_invariants.npz`: full verified response vectors, invariants, critical-axis arrays, and 721-angle terminal observer scans.
- `results/channel_vector_angle_and_dot.png`: mutual angle and normalized dot.
- `results/channel_oriented_area_and_handedness.png`: signed area and normalized orientation.
- `results/observer_critical_angles_and_compensation.png`: observer-zero axes and observer-dependent cancellation.
- `results/report.json` and `results/SCIENTIFIC_REPORT_SERVER.md`: numerical metrics, checks and limitations.
- `SHA256SUMS_OUTPUT.json`: executable output hashes.

## Reproducibility and limitations

See `MATHEMATICAL_PROTOCOL.md` and `SCIENTIFIC_DECISION_RULES.md`. Signed area is SO(2), but not O(2), invariant; source-channel attribution is conditional on model coordinates. A zero projected component is not a physical zero vector or a return event. All results are finite-interval numerical diagnostics, not certified material memory. The seed-invariant geometrical trends are independent of whether the same physical events are labeled by xi or tau, but do not give a new fixed-tau parameter derivative.
