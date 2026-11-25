#!/usr/bin/env bash
set -euo pipefail
BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PYTHON_BIN:-python3.13}"
cd "$BUNDLE"
"$PY" scripts/verify_manifest.py
"$PY" tests/test_factorial_geometry.py
"$PY" scripts/metric_path_attribution.py --preflight
case "${1:-}" in
    --preflight-only) exit 0 ;;
    --smoke-only) OPENBLAS_NUM_THREADS=1 "$PY" scripts/metric_path_attribution.py --smoke;exit 0;;
    ''|--require-parent|--full) ;;
    *) echo "Unknown argument $1" >&2;exit 2;;
esac
OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}" "$PY" scripts/metric_path_attribution.py --full
"$PY" scripts/package_sendback.py
