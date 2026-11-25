#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
bash scripts/run_on_server.sh "$@"
