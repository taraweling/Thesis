#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

JOB_SCRIPT="run_walker.sbatch"

if ! command -v sbatch >/dev/null 2>&1; then
  echo "Error: sbatch not found. Run this on the HPC login node." >&2
  exit 1
fi

if [[ ! -f "$JOB_SCRIPT" ]]; then
  echo "Error: ${JOB_SCRIPT} not found in ${SCRIPT_DIR}" >&2
  exit 1
fi

mkdir -p logs

submit_output="$(sbatch "$JOB_SCRIPT")"
echo "$submit_output"

job_id="$(echo "$submit_output" | awk '{print $4}')"
if [[ -n "${job_id:-}" ]]; then
  echo "Queued job ID: ${job_id}"
  squeue -j "$job_id" || true
fi
