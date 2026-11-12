#!/usr/bin/env bash
#SBATCH --job-name=aqua2-analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=aqua2-%j.log
set -euo pipefail
: "${MATLAB_BIN:?Set the MATLAB executable path}"
: "${AQUA2_RUNNER:?Set the externally installed batch-runner function name}"
: "${AQUA2_INPUT:?Set an input TIFF path}"
: "${AQUA2_PRESET:?Set an explicit preset CSV path}"
# Submit from this repository root. No account, partition or institution is assumed.
"$MATLAB_BIN" -batch "setup_aqua2_tools; run_aqua2_from_environment"
