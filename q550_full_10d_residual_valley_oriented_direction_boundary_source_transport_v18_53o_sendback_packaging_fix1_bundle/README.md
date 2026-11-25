# Q550 v18.53o — Sendback packaging recovery fix1

## Scope

The physical audit completed all **10 IVPs**, including four tighter integrations, and reported:

`AUTHENTIC_ORIENTED_10D_VALLEY_AND_BOUNDARY_GEOMETRY__NO_MATERIAL_HOLONOMY`

The original `scripts/package_sendback.py` then failed because it combined a `str` and `PosixPath` at its `z.write` expression. It had **already written the `OUTPUT_MANIFEST.json`** before the Python exception.

This recovery is **packaging-only**. It does not run the physical solver, edit the original bundle, regenerate science, or change the existing source and output manifest contents. It creates a verified sendback ZIP outside the original bundle directory using `Path(root.name) / p.relative_to(root)` for archive member names.

## How to run on the existing server results

After extracting this fix1 patch alongside the original v18.53o bundle in `/root/elab1/hysteresis`:

```bash
cd /root/elab1/hysteresis
PATCH=q550_full_10d_residual_valley_oriented_direction_boundary_source_transport_v18_53o_sendback_packaging_fix1_bundle
unzip -q "${PATCH}.zip"
PYTHON_BIN=python3.13 bash "$PATCH/run_recovery.sh"
```

To use a different server workspace, pass its parent directory to `run_recovery.sh`.

The recovery overwrites only the previously failed **sendback ZIP**, leaving original physical files, source scripts, result data and original manifests untouched. For safety, the Python helper itself requires an explicit `--overwrite` when an output file already exists; the shell wrapper opts into this because the failed packager may have created an incomplete ZIP.

Expected last lines:

```
PHYSICAL_RESULTS_PRESENT 10 IVPS; 4 TIGHT; 846+846 CURVE ROWS
SOURCE_MANIFEST_PASS 14
OUTPUT_MANIFEST_PASS <count>
ZIP_CRC_AND_BOTH_MANIFESTS_PASS <count> members
SEND_BACK_RECOVERED /root/elab1/hysteresis/q550_full_10d_residual_valley_oriented_direction_boundary_source_transport_v18_53o_bundle_server_sendback.zip
SHA256 <digest>
NO_PHYSICAL_IVPS_EXECUTED
```

Send that `..._server_sendback.zip` to ChatGPT for scientific review.

## Validation and scope limits

- The recovery checks the successful scientific decision and 10/4 IVP counts in `results/report.json`.
- It verifies expected fixed/valley rows (846 each), six per-case records and existing NPZ/CSV/report artifacts.
- It independently verifies **both** `SOURCE_MANIFEST.json` and `OUTPUT_MANIFEST.json` against on-disk files.
- It writes a temporary ZIP, checks every member's CRC, checks all referenced hashes *inside the ZIP*, then atomically replaces the failed sendback archive.
- It does **not** certify the physical conclusions. Only the original server experiment can provide those results.
- The old unpatched `scripts/package_sendback.py` remains a historical source file inside the recovered sendback; this small independent recovery script is the documented patch, not a rewrite of the original genealogy.
- On a missing or incomplete physical result, the helper **fails closed** instead of fabricating an output manifest or re-running the solver.

The complete input v18.53o distribution ZIP is not needed to execute this recovery; the original extracted directory with successful results is required.
