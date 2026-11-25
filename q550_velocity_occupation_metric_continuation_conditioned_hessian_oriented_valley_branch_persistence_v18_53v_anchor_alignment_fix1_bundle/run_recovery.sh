#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
PATCH=q550_velocity_occupation_metric_continuation_conditioned_hessian_oriented_valley_branch_persistence_v18_53v_anchor_alignment_fix1_bundle
export PYTHON_BIN="${PYTHON_BIN:-python3.13}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
"$PYTHON_BIN" "$PATCH/recover_v18_53v.py" --self-test
"$PYTHON_BIN" "$PATCH/recover_v18_53v.py" --recover
