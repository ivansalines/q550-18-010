#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
case "${1:---require-parent}" in
  --require-parent|--full) MODE=full ;;
  --preflight-only) MODE=preflight ;;
  --smoke) MODE=smoke ;;
  *) echo "Unsupported argument $1" >&2; exit 2 ;;
esac
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
"$PYTHON_BIN" scripts/verify_manifest.py
"$PYTHON_BIN" tests/test_continuation_geometry.py
"$PYTHON_BIN" tests/test_packaging.py
"$PYTHON_BIN" scripts/metric_continuation.py --preflight
if [[ "$MODE" == "preflight" ]]; then
  echo 'V18_53V_PREFLIGHT_COMPLETE__NO_PHYSICAL_IVPS'
  exit 0
fi
if [[ "$MODE" == "smoke" ]]; then
  "$PYTHON_BIN" scripts/metric_continuation.py --smoke
  exit 0
fi
"$PYTHON_BIN" scripts/metric_continuation.py --full
"$PYTHON_BIN" scripts/package_sendback.py
