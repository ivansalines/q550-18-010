# Q550 — v18.53v Anchor-Alignment Fix1 & Sendback Recovery

The 204-node metric-continuation calculation completed, but the original scientific run stopped at the final `PARENT_U_REGRESSION_BREACH 0.000328778` condition. That final gate incorrectly compared arclengths of *different anchored curves*:

- Historical v18.53u observed its 141-point curve around the **v18.53t historical curve center**.
- v18.53v observed its 141-point curve around the **optimized canonical minimum from v18.53u**.

For example, original/shared: historical curve anchor = `14.328678903551179`; new v curve anchor = `14.327390213582948`. This alters arclength although minimizer coordinates and residual norms agree to numerical precision. The difference must not be mistaken for a branch failure or corrected by merely increasing the threshold.

## Recover on server, without rerunning the 204-node post-processing

```bash
cd /root/elab1/hysteresis
P=q550_velocity_occupation_metric_continuation_conditioned_hessian_oriented_valley_branch_persistence_v18_53v_anchor_alignment_fix1_bundle
unzip -q "$P.zip"
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 bash "$P/run_recovery.sh" \
  2>&1 | tee v18_53v_anchor_alignment_fix1_recovery.log
```

The original v18.53v extracted folder must still contain the 204-node result files that were saved before the failed gate. If they do not exist, the script fails explicitly rather than inventing missing results.

The recovery:

1. Verifies the original source manifest (14/14), immutable v18.53u and v18.53t nested ZIP parent manifests and their hashes.
2. Checks **204 metric-continuation rows**, **28,764 curve rows**, **36 raw parent checks**, **12 family rows**, and **204 oriented 45D polygons**.
3. Uses the authenticated v18.53u width-search algorithm (25-point bracket search) on the *same 141-point origin interval* as v18.53v (whose original solver uses 21-point brackets). At matched anchors, all historical min fields and arclengths must satisfy explicit tolerances.
4. Preserves the original `metric_lambda_parent_u_replay.csv` to document the incorrectly compared paths, and adds `metric_lambda_parent_u_matched_anchor.csv`.
5. Generates final `report.json` and `SCIENTIFIC_REPORT_SERVER.md` based on already computed node diagnostics, then invokes the original v18.53v ZIP packager to make its sendback with CRC and output SHA-256 checks.

**No original source code, archived parent ZIP, physical ODE, or original CSV is modified.** The added files are a provenance-aware recovery report, not a new scientific experiment. The output ZIP has the original full name:

`q550_velocity_occupation_metric_continuation_conditioned_hessian_oriented_valley_branch_persistence_v18_53v_bundle_server_sendback.zip`

## Scientific limits

This is conditioned *numerical metric-family continuation*, not a certificate of continuum flow, material holonomy, or full-state 11D recurrence. Metric-weighted 45D orientation comparisons use a specific SPD square-root identification. The far endpoints remain in formal-clock continuation beyond xi=25.6.
