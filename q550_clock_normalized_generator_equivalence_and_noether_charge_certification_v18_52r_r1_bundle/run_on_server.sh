#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_manifest.py
"$PYTHON_BIN" tests/test_synthetic.py
# Exact symbolic proof runs with SymPy when available; Fraction fallback needs no SymPy.
"$PYTHON_BIN" scripts/symbolic_noether.py
set +e
"$PYTHON_BIN" scripts/clock_audit.py --strict "$@"
STATUS=$?
set -e
"$PYTHON_BIN" scripts/package_sendback.py
exit "$STATUS"
