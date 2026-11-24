#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_source_manifest.py
"$PYTHON_BIN" tests/test_synthetic.py
"$PYTHON_BIN" scripts/exact_symbolic_proof.py
set +e
"$PYTHON_BIN" scripts/run_audit.py --mode server "$@"
AUDIT_RC=$?
set -e
"$PYTHON_BIN" scripts/pack_sendback.py
exit "$AUDIT_RC"
