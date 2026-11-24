# Q550 v18.51y — Scientific Report (interpretive note)

The source-support and coordinate-covariance diagnostics are replayable from the authenticated v18.51t/v18.51w/v18.51x inputs. The original channel decomposition is unique **conditional on fixed historical state coordinates, fixed channel labels, and disjoint source supports**. It is not a globally unique physical actuator decomposition.

The correctly transformed tangent generator and source projectors preserve channel identity under passive time-dependent chart changes. Omission of projector transport or the frame connection creates a purely representational discrepancy. The off-support Jacobian column ratios measure coordinate-dependent instantaneous coupling, not invariant energy.

The current numerical outcome, scientific check ledger, and source/output hashes are generated in `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, and `SHA256SUMS_OUTPUT.json` by the runner. Independent server replay remains the reproducibility check; no nonlinear model IVPs are performed at this stage.
