#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_manifest.py || exit 4
"$PYTHON_BIN" tests/test_synthetic.py || exit 5
"$PYTHON_BIN" scripts/run_audit.py "$@"
status=$?
"$PYTHON_BIN" scripts/package_sendback.py || exit 6
exit "$status"
