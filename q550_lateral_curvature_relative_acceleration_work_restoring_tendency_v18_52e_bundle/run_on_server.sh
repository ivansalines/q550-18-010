#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON_BIN:-python3.13}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
"$PY" scripts/verify_source_manifest.py
"$PY" tests/test_synthetic.py
"$PY" scripts/run_audit.py --mode server "$@"
"$PY" tests/test_results.py
"$PY" scripts/write_output_manifest.py
"$PY" scripts/pack_sendback.py
