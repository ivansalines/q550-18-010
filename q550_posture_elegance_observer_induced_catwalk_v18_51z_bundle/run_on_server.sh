#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
P="${PYTHON_BIN:-python3.13}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
"$P" scripts/verify_source_manifest.py
"$P" tests/test_synthetic.py
"$P" scripts/run_audit.py
"$P" tests/test_results.py
"$P" scripts/write_output_manifest.py
"$P" scripts/pack_sendback.py
