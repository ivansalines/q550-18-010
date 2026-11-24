#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_manifest.py
"$PYTHON_BIN" tests/test_synthetic.py
"$PYTHON_BIN" tests/test_parameter_lineage.py
"$PYTHON_BIN" scripts/symbolic_energy.py
set +e
"$PYTHON_BIN" scripts/run_audit.py "$@"
STATUS=$?
set -e
"$PYTHON_BIN" scripts/package_sendback.py
exit "$STATUS"
