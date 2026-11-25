#!/usr/bin/env bash
set -euo pipefail
BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PYTHON_BIN:-python3.13}"
cd "$BUNDLE"
mkdir -p results
"$PY" scripts/verify_manifest.py
"$PY" tests/test_synthetic.py
"$PY" tests/test_packaging.py
OPENBLAS_NUM_THREADS=1 "$PY" scripts/run_audit.py --preflight
case "${1:---require-parent}" in
  --preflight-only) exit 0 ;;
  --smoke-only) OPENBLAS_NUM_THREADS=1 "$PY" scripts/run_audit.py --smoke; exit 0 ;;
  --require-parent|--full) ;;
  *) echo "Unknown argument: $1" >&2; exit 2 ;;
esac
OPENBLAS_NUM_THREADS=1 "$PY" scripts/run_audit.py --full
"$PY" scripts/package_sendback.py
