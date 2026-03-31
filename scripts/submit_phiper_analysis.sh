#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./submit_phiper_analysis.sh <script> <project> <group> [all] [longitudinal] [parquet_name] [max_mem] [cpus] [time_limit]

Arguments:
  script           Path to the R analysis script
  project          Project directory, e.g. IBD-Berlin
  group            Active group name, e.g. group_test

Optional arguments:
  all              TRUE/FALSE, whether to run the "all" analysis block
                   default: FALSE

  longitudinal     TRUE/FALSE, default longitudinal mode
                   default: FALSE

  parquet_name     Name of parquet file inside <project>/Data/
                   default: <project>.parquet

  max_mem          Slurm memory request
                   default: 40G

  cpus             Slurm CPUs per task
                   default: 30

  time_limit       Slurm walltime
                   default: 05:00:00

Examples:
  ./submit_phiper_analysis.sh src/R/02-run_phiper_analysis.R IBD-Berlin group_test

  ./submit_phiper_analysis.sh src/R/02-run_phiper_analysis.R IBD-Berlin group_test TRUE FALSE

  ./submit_phiper_analysis.sh src/R/02-run_phiper_analysis.R IBD-Berlin group_test TRUE FALSE IBD-Berlin_newlabel.parquet

  ./submit_phiper_analysis.sh src/R/02-run_phiper_analysis.R IBD-Berlin pre-screening TRUE FALSE IBD-Berlin.parquet 60G 16 24:00:00

Notes:
  - The job is submitted with sbatch.
  - Log files go to logs/<project>_<group>_<jobid>.out/.err
  - parquet_name is interpreted relative to <project>/Data/
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 3 ]]; then
  usage
  exit 1
fi

script="${1:?Need R script}"
project="${2:?Need project dir}"
group="${3:?Need active group}"
all="${4:-FALSE}"
longitudinal="${5:-FALSE}"
parquet_name="${6:-${project}.parquet}"
max_mem="${7:-40G}"
cpus="${8:-30}"
time_limit="${9:-05:00:00}"

log_file="${project}_${group}"

mkdir -p logs

echo "Submitting job with:"
echo "  script        = ${script}"
echo "  project       = ${project}"
echo "  group         = ${group}"
echo "  all           = ${all}"
echo "  longitudinal  = ${longitudinal}"
echo "  parquet_name  = ${parquet_name}"
echo "  max_mem       = ${max_mem}"
echo "  cpus          = ${cpus}"
echo "  time_limit    = ${time_limit}"

sbatch \
  --job-name="${project}" \
  --output="logs/${log_file}_%j.out" \
  --error="logs/${log_file}_%j.err" \
  --mem="${max_mem}" \
  --cpus-per-task="${cpus}" \
  --time="${time_limit}" \
  scripts/run_phiper.sh \
  "${script}" \
  "${project}" \
  "${group}" \
  "${all}" \
  "${longitudinal}" \
  "${parquet_name}" 
