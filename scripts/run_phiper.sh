#!/bin/bash
#SBATCH --job-name=phiper
#SBATCH --output=logs/IBD-Chile_%x_%j.out
#SBATCH --error=logs/IBD-Chile_%x_%j.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=40G

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE:-unknown} MB"

module purge
module load R
module load Pandoc/3.9

SCRIPT="${1:?Usage: sbatch script.sh <R_SCRIPT> <PROJECT_DIR> <ACTIVE_GROUP> [ALL] [DEFAULT_LONGITUDINAL] [PARQUET_NAME]}"
PROJECT_DIR="${2:?Need PROJECT_DIR}"
ACTIVE_GROUP="${3:?Need ACTIVE_GROUP}"
ALL_FLAG="${4:-FALSE}"
DEFAULT_LONGITUDINAL="${5:-FALSE}"
PARQUET_NAME="${6:-${PROJECT_DIR}.parquet}"

if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  MAX_GB_ARG=$(( SLURM_MEM_PER_NODE / 1024 - 1 ))
  if [[ "${MAX_GB_ARG}" -lt 1 ]]; then
    MAX_GB_ARG=1
  fi
else
  MAX_GB_ARG=40
fi

echo "Using script: ${SCRIPT}"
echo "Using PROJECT_DIR: ${PROJECT_DIR}"
echo "Using ACTIVE_GROUP: ${ACTIVE_GROUP}"
echo "Using ALL: ${ALL_FLAG}"
echo "Using DEFAULT_LONGITUDINAL: ${DEFAULT_LONGITUDINAL}"
echo "Using PARQUET_NAME: ${PARQUET_NAME}"
echo "Using MAX_GB: ${MAX_GB_ARG}"

Rscript --no-save --no-restore "${SCRIPT}" \
  N_CORES="${SLURM_CPUS_PER_TASK}" \
  MAX_GB="${MAX_GB_ARG}" \
  LOG=true \
  FORCE=false \
  ACTIVE_GROUP="${ACTIVE_GROUP}" \
  ALL="${ALL_FLAG}" \
  DEFAULT_LONGITUDINAL="${DEFAULT_LONGITUDINAL}" \
  PROJECT_DIR="${PROJECT_DIR}" \
  PARQUET_NAME="${PARQUET_NAME}"
