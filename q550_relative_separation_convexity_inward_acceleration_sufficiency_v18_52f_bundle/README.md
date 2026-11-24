# Q550 — v18.52f: Relative-Separation Convexity & Inward-Acceleration Sufficiency Audit

This is a reproducible, **nonintegrating** audit of eight genuine finite-kick physical trajectories from the authenticated v18.52e server sendback. The complete original parent server archive is embedded in `inputs/`; it includes the original field provenance and prior run scripts.

The v18.52f script independently rebuilds the three signed scalar kinematic diagnostics from archived full five-dimensional physical Cartesian vectors at the same elapsed material time:

- `positive = ||delta V||^2`;
- `inward = -delta Q · delta a`;
- `convexity = positive - inward = d²[||delta Q||²/2]/d tau²`.

It evaluates the necessary **convexity reversal threshold** `inward/positive >= 1` separately from the **separation-rate reversal threshold** `G=delta Q·delta V <= 0`. A scalar counterfactual amplification factor `positive/inward` quantifies the instantaneous shortfall; it is **not** a physical parameter, modified force, or dynamical continuation.

## Server run

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_relative_separation_convexity_inward_acceleration_sufficiency_v18_52f_bundle.zip
cd q550_relative_separation_convexity_inward_acceleration_sufficiency_v18_52f_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52f_server.log
```

No expensive nonlinear integration is required: all genuine trajectories were integrated on the user's server in v18.52e. The new audit is a distinct, nontrivial kinematic sufficiency test with independent spline and quadrature checks.

The runner verifies authentic source hashes and CRC, performs synthetic adversarial tests, processes the original vectors, creates CSV, figure, JSON/scientific reports, checks outputs, writes the SHA-256 manifest and packs the server sendback.

`results/report.json` reports **sampled-domain** properties only. Positive sample values or shape-preserving scalar interpolants do NOT constitute rigorous continuous bounds for the nonlinear flow; neither a loss of convexity nor a positional return is certified outside the observed horizon.

Requirements: Python >=3.10, NumPy, SciPy and Matplotlib. Text and code comments are in English. Full source genealogy stays inside the parent sendback.
