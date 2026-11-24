#!/usr/bin/env bash
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts:$ROOT/reference:${PYTHONPATH:-}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
PY="${PYTHON_BIN:-python3.13}"
"$PY" scripts/verify_manifest.py || exit 11
"$PY" tests/test_synthetic.py || exit 12
"$PY" tests/test_archived_rhs_crosscheck.py || exit 13
"$PY" scripts/run_audit.py "$@"
rc=$?
"$PY" scripts/package_sendback.py
exit "$rc"
