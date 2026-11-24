#!/usr/bin/env bash
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts:$ROOT/reference:${PYTHONPATH:-}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
PY="${PYTHON_BIN:-python3.13}"
"$PY" scripts/verify_manifest.py
rc=$?
if [[ $rc -eq 0 ]]; then "$PY" tests/test_synthetic.py; rc=$?; fi
if [[ $rc -eq 0 ]]; then "$PY" tests/test_archived_rhs_crosscheck.py; rc=$?; fi
if [[ $rc -eq 0 ]]; then "$PY" scripts/run_audit.py "$@"; rc=$?; fi
"$PY" scripts/package_sendback.py
exit "$rc"
