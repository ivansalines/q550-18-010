#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_manifest.py
"$PYTHON_BIN" tests/test_synthetic.py
"$PYTHON_BIN" scripts/run_audit.py --preflight
"$PYTHON_BIN" scripts/run_audit.py "$@"
"$PYTHON_BIN" scripts/package_sendback.py
