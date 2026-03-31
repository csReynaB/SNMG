#!/bin/bash
#SBATCH --job-name=IBD-report
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

module purge
module load R
module load Quarto
module load Pandoc/3.9

SCRIPT="src/R/03-render_phiper_reports.R"

BASE_DIR_REL="${1:-}"
GROUP_COLS="${2:-}"

if [[ -z "$BASE_DIR_REL" ]]; then
  echo "Usage: sbatch scripts/get_reports.sh <relative_or_absolute_base_dir> <group1,group2,...>"
  exit 1
fi

if [[ -z "$GROUP_COLS" ]]; then
  echo "Usage: sbatch scripts/get_reports.sh <relative_or_absolute_base_dir> <group1,group2,...>"
  exit 1
fi

BASE_DIR_ABS="$(realpath "$BASE_DIR_REL")"

echo "Working directory: $(pwd)"
echo "Base dir: $BASE_DIR_ABS"
echo "Group cols: $GROUP_COLS"

Rscript --no-save --no-restore "$SCRIPT" \
  BASE_DIR="$BASE_DIR_ABS" \
  GROUP_COLS="$GROUP_COLS"
