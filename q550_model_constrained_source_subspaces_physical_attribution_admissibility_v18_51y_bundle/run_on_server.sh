#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
"$PYTHON_BIN" scripts/verify_source_manifest.py
"$PYTHON_BIN" tests/test_model_constrained_synthetic.py
"$PYTHON_BIN" scripts/run_audit.py
"$PYTHON_BIN" tests/test_results.py
"$PYTHON_BIN" scripts/write_output_manifest.py
"$PYTHON_BIN" scripts/pack_sendback.py
