#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
: "${PYTHON_BIN:=python3.13}"
: "${OPENBLAS_NUM_THREADS:=1}"
export OPENBLAS_NUM_THREADS
export OMP_NUM_THREADS=1
"$PYTHON_BIN" scripts/verify_source_manifest.py
"$PYTHON_BIN" tests/test_synthetic.py
"$PYTHON_BIN" tests/test_mock_clock_ivp.py
"$PYTHON_BIN" scripts/run_audit.py --mode preflight
"$PYTHON_BIN" scripts/run_audit.py --mode server "$@"
"$PYTHON_BIN" tests/test_results.py
"$PYTHON_BIN" scripts/write_output_manifest.py
"$PYTHON_BIN" scripts/pack_sendback.py
