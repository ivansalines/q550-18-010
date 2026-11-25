#!/usr/bin/env bash
set -euo pipefail
PARENT_DIR="${1:-/root/elab1/hysteresis}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B=q550_full_10d_residual_valley_oriented_direction_boundary_source_transport_v18_53o_bundle
PYTHON_BIN="${PYTHON_BIN:-python3.13}"
exec "$PYTHON_BIN" "$HERE/recover_q550_v18_53o_sendback.py" \
  --bundle-dir "$PARENT_DIR/$B" --overwrite
