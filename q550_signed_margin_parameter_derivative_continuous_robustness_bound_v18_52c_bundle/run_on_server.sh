#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
"$PYTHON_BIN" scripts/verify_source_manifest.py
"$PYTHON_BIN" tests/test_synthetic.py
status=0
"$PYTHON_BIN" scripts/run_audit.py "$@" || status=$?
if [[ -f results/report.json ]]; then
    "$PYTHON_BIN" tests/test_results.py || status=4
fi
"$PYTHON_BIN" scripts/write_output_manifest.py
"$PYTHON_BIN" scripts/pack_sendback.py
if [[ "$status" -ne 0 ]]; then
    echo "Q550 v18.52c: audit partial or failed; sendback was still packed; exit=$status" >&2
fi
exit "$status"
