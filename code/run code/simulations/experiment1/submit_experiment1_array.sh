#!/bin/bash
#SBATCH --job-name=experiment1
#SBATCH --array=1-200
#SBATCH --time=00:40:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-experiment1-%A_%a.out

set -eo pipefail

# Required: absolute project path (the directory containing code/ and output/).
# Submit with: sbatch "/absolute/path/to/project/code/run code/simulations/experiment1/submit_experiment1_array.sh"
# Supply account/partition/QOS as sbatch options only if your cluster requires them.
PROJECT_ROOT=""

# Leave blank if Rscript and matlab are already available in the job environment.
# Otherwise, fill in the module names provided by your cluster.
R_MODULE=""
MATLAB_MODULE=""

if [[ -z "${PROJECT_ROOT}" || "${PROJECT_ROOT}" != /* ]]; then
  echo "Set PROJECT_ROOT to the absolute project directory in this submit script." >&2
  exit 1
fi

export TREEFA_PROJECT_ROOT="${PROJECT_ROOT}"
export TREEFA_SIM_DIR="${PROJECT_ROOT}/code/run code/simulations/experiment1"
export TREEFA_CLUSTER_CODE_DIR="${PROJECT_ROOT}/code/core_code"
export TREEFA_R_LIB="${TREEFA_CLUSTER_CODE_DIR}/rlib"
export TREEFA_EXPERIMENT_NAME=experiment1
export TREEFA_RESULT_DIR="${PROJECT_ROOT}/output/${TREEFA_EXPERIMENT_NAME}"
RUN_SCRIPT="${TREEFA_SIM_DIR}/run_cluster_experiment1.R"

for required_dir in "${TREEFA_SIM_DIR}" "${TREEFA_CLUSTER_CODE_DIR}"; do
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
# runtime output to the original per-task log filenames under output/experiment1/.
exec > "${TREEFA_RESULT_DIR}/results_${SLURM_ARRAY_TASK_ID}.Rout" \
     2> "${TREEFA_RESULT_DIR}/errors_${SLURM_ARRAY_TASK_ID}.err"

for module_name in "${R_MODULE}" "${MATLAB_MODULE}"; do
  if [[ -n "${module_name}" ]]; then
    if ! command -v module >/dev/null 2>&1; then
      echo "The module command is unavailable. Make Rscript and matlab available before submitting." >&2
      exit 1
    fi
    module load "${module_name}"
  fi
done

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found. Load R before submitting or fill in R_MODULE." >&2
  exit 1
fi

if [[ -z "${MATLAB_BIN:-}" ]]; then
  MATLAB_BIN="$(command -v matlab || true)"
fi
if [[ ! -f "${MATLAB_BIN}" || ! -x "${MATLAB_BIN}" ]]; then
  echo "MATLAB not found. Load MATLAB, fill in MATLAB_MODULE, or export MATLAB_BIN with its executable path." >&2
  exit 1
fi
export MATLAB_BIN

export TREE_VARIANT_RS_SELECTION=validation
export TREE_VARIANT_RS_ADAPTIVE_METHODS=
export TREE_VARIANT_RS_GAMMA_NLAM=50
export TREE_VARIANT_RS_GAMMA_MIN_RATIO=1e-4
export TREEFA_KEEP_RS_FILES=false

export TREEFA_RS_WORKDIR="${TREEFA_RESULT_DIR}/rs_runs/array_${SLURM_ARRAY_TASK_ID}"

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "Using Rscript: $(command -v Rscript)"
echo "Using MATLAB_BIN=${MATLAB_BIN}"
echo "TREEFA_SIM_DIR=${TREEFA_SIM_DIR}"
echo "TREEFA_CLUSTER_CODE_DIR=${TREEFA_CLUSTER_CODE_DIR}"
echo "TREEFA_RESULT_DIR=${TREEFA_RESULT_DIR}"

cd "${TREEFA_SIM_DIR}"
exec Rscript --vanilla "${RUN_SCRIPT}"
