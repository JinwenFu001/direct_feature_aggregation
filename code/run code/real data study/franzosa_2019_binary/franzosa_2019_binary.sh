#!/bin/bash
#SBATCH --output=franzosa_2019_binary_%A_%a.slurm.out
#SBATCH --array=1-200
#SBATCH --job-name=franzosa_2019_binary
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=06:00:00

set -eo pipefail

# Required: absolute project path (the directory containing code/, data/ and output/).
# Submit with: sbatch "/absolute/path/to/project/code/run code/real data study/franzosa_2019_binary/franzosa_2019_binary.sh"
# Supply account/partition/QOS as sbatch options only if your cluster requires them.
PROJECT_ROOT=""

# Leave blank if Rscript is already available in the job environment.
# Otherwise, fill in the R module name provided by your cluster.
R_MODULE=""

if [[ -z "${PROJECT_ROOT}" || "${PROJECT_ROOT}" != /* ]]; then
  echo "Set PROJECT_ROOT to the absolute project directory in this submit script." >&2
  exit 1
fi

export TREEFA_PROJECT_ROOT="${PROJECT_ROOT}"
export TREEFA_STUDY_DIR="${PROJECT_ROOT}/code/run code/real data study/franzosa_2019_binary"
export TREEFA_CLUSTER_CODE_DIR="${PROJECT_ROOT}/code/core_code"
export TREEFA_R_LIB="${TREEFA_CLUSTER_CODE_DIR}/rlib"
export TREEFA_REAL_DATA_DIR="${PROJECT_ROOT}/data"
export TREEFA_RESULT_DIR="${PROJECT_ROOT}/output/real data study/franzosa_2019_binary"
RUN_SCRIPT="${TREEFA_STUDY_DIR}/franzosa_2019_binary.R"

for required_dir in "${TREEFA_STUDY_DIR}" "${TREEFA_CLUSTER_CODE_DIR}" "${TREEFA_REAL_DATA_DIR}"; do
  if [[ ! -d "${required_dir}" ]]; then
    echo "Required directory not found: ${required_dir}" >&2
    exit 1
  fi
done
if [[ ! -r "${RUN_SCRIPT}" ]]; then
  echo "R script not found or not readable: ${RUN_SCRIPT}" >&2
  exit 1
fi
if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
  echo "Submit this script with sbatch; SLURM_ARRAY_TASK_ID must be a positive integer." >&2
  exit 1
fi

mkdir -p "${TREEFA_RESULT_DIR}"

# SBATCH paths cannot expand PROJECT_ROOT, and Slurm opens its log before this
# script runs. Keep the bootstrap log in the submission directory, then send
# combined runtime output to the original per-task Rout filename under
# output/real data study/franzosa_2019_binary/.
exec > "${TREEFA_RESULT_DIR}/franzosa_2019_binary_${SLURM_ARRAY_TASK_ID}.Rout" 2>&1

if [[ -n "${R_MODULE}" ]]; then
  if ! command -v module >/dev/null 2>&1; then
    echo "The module command is unavailable. Make Rscript available before submitting." >&2
    exit 1
  fi
  module load "${R_MODULE}"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found. Load R before submitting or fill in R_MODULE." >&2
  exit 1
fi

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "Using Rscript: $(command -v Rscript)"
echo "TREEFA_STUDY_DIR=${TREEFA_STUDY_DIR}"
echo "TREEFA_CLUSTER_CODE_DIR=${TREEFA_CLUSTER_CODE_DIR}"
echo "TREEFA_R_LIB=${TREEFA_R_LIB}"
echo "TREEFA_REAL_DATA_DIR=${TREEFA_REAL_DATA_DIR}"
echo "TREEFA_RESULT_DIR=${TREEFA_RESULT_DIR}"

cd "${TREEFA_STUDY_DIR}"
exec Rscript --vanilla "${RUN_SCRIPT}"
