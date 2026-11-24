#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
: "${PYTHON_BIN:=python3.13}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
"$PYTHON_BIN" scripts/verify_lineage.py
"$PYTHON_BIN" tests/test_synthetic.py
"$PYTHON_BIN" scripts/run_audit.py --mode preflight "$@"
"$PYTHON_BIN" scripts/run_audit.py --mode audit "$@"
"$PYTHON_BIN" scripts/write_manifest.py
"$PYTHON_BIN" scripts/pack_sendback.py
