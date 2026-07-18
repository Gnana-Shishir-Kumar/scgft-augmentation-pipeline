#!/usr/bin/env bash
# Run the scGFT augmentation pipeline stages in order. Stops on first failure.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found on PATH. Install R and ensure Rscript is available." >&2
  exit 1
fi

run_stage() {
  local stage_name="$1"
  local script_path="$2"

  echo ""
  echo "============================================================"
  echo "STAGE: ${stage_name}"
  echo "SCRIPT: ${script_path}"
  echo "============================================================"

  if [[ ! -f "${script_path}" ]]; then
    echo "ERROR: Stage failed: ${stage_name} — missing script ${script_path}" >&2
    exit 1
  fi

  if ! Rscript "${script_path}"; then
    echo "" >&2
    echo "ERROR: Stage failed: ${stage_name} (${script_path})" >&2
    echo "Pipeline stopped. Fix the error above and re-run." >&2
    exit 1
  fi

  echo "STAGE OK: ${stage_name}"
}

run_stage "01_download_data" "R/01_download_data.R"
run_stage "02_preprocess"    "R/02_preprocess.R"
run_stage "03_run_scgft"     "R/03_run_scgft.R"
run_stage "04_export_json"   "R/04_export_json.R"

echo ""
echo "============================================================"
echo "Pipeline finished successfully."
echo "JSON outputs are under: ${ROOT}/results/"
echo "============================================================"
