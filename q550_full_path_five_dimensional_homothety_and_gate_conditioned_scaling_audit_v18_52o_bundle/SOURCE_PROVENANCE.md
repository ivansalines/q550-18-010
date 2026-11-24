# Source provenance — v18.52o

## Verified physical source

`v18.52l/results/extended_horizon_5d_vectors_and_traces.npz` is pinned by SHA-256 `d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`. It contains eight physical 5D relative position/velocity trajectories, physical material time and reference `D/G` calculations. The bundle does not embed or recompute this heavy ancestor.

## Observer ancestry

`reference/PARENT_v18_52m_report.json` stores the historical 5D PCA direction and the already-audited first-half training convention. `reference/PARENT_v18_52n_report.json` stores the exact v18.52n interpretation and matched gate slope. Historical `m≈-0.308` is an observational hypothesis because the user has not supplied the four exact 2D point coordinates. It must not be represented as a measured best-fit from the v18.52l section roots (those were `Q2=0`). The angle `atan(-0.308)≈-17.119°` is distinct from the deliberately chosen observer tilt `±17.04°`.

## Reproducibility contract

Actual generator: `scripts/run_audit.py`. Entry script: `run_on_server.sh`. Synthetic methods tests: `tests/test_synthetic.py`. Historical verification: `scripts/verify_lineage.py`. Scripted SHA-256/ZIP packing: `scripts/write_manifest.py`, `scripts/pack_sendback.py`. `requirements.txt` lists dependencies; source code and SHA-256 manifests preserve the executable scientific gesture, not merely derived plots or CSVs. The sendback contains all of these.

## Caveat

No scientific result is transferred from the local development smoke; all real numbers must be read from the forthcoming full-resolution server sendback.
