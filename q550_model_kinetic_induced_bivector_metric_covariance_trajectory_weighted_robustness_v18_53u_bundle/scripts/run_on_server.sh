#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
case "${1:-}" in
  --require-parent|"") ;;
  *) echo "Unsupported flag: $1" >&2; exit 2 ;;
esac
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_manifest.py
"$PYTHON_BIN" tests/test_metric_geometry.py
"$PYTHON_BIN" scripts/run_audit.py --preflight
"$PYTHON_BIN" scripts/run_audit.py --full
"$PYTHON_BIN" scripts/package_sendback.py
