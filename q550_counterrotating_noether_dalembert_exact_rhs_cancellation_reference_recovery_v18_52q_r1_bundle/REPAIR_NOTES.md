# v18.52q-r1 — missing external source reference repair

The original v18.52q runner required an external **SHA-pinned v18.51s server
sendback ZIP** and stopped the run whenever that archive was absent. The
original v18.52q executable and report were not a self-contained current-engine
certificate. The parent `v18.52l` 10D NPZ is also an external reference.

## What this revision changes

- `--mode auto` is the default. Runs exact historic symbolic/counterfactual/RHS
  probes regardless of whether the v18.51s sendback is available.
- Native 10D traces are verified only if the genuine v18.52l file is present;
  its SHA-256 is mandatory.
- The current-engine `structural_jet` probes run **only** if an authentic
  v18.51s archive is present. Its SHA-256 must match; no fabricated fallback.
- Missing sources are explicitly listed in `results/report.json` and console
  output. The certificate status is reduced, not silently upgraded.
- `--mode server` remains a strict option and returns nonzero for incomplete
  current-engine certification. `--mode archived` explicitly omits current engine.

## How to restore the complete independent check

Place `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`
in `/root/elab1/hysteresis`, or supply
`bash run_on_server.sh --source-s-path /absolute/path/to/source.zip --mode server`.
The expected SHA-256 is recorded in `inputs_config.json`.

This revision is not a new dynamical model and introduces no new integrations.
The algebraic cancellation in the historical documented equations remains exact.
