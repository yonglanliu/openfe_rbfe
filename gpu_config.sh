#!/bin/bash

#SBATCH --job-name=rbfe
#SBATCH --partition=gpu 
#SBATCH --time=120:00:00
#SBATCH --gres=gpu:v100x:3
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --array=0-17
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err    

# -------- User settings --------
OPENFE_ENV="openfe"               # conda env name
JSON_DIR="/vf/users/liuy48/fep/openfe/TYK2/json"     # folder containing 18 *.json files
OUT_ROOT="/vf/users/liuy48/fep/openfe/TYK2/results"        # where replicate_* folders will be created
N_REP=1                                # replicates per transformation
# --------------------------------

mkdir -p logs
for ((i=0; i<${N_REP}; i++)); do
    mkdir -p "${OUT_ROOT}/replicate_${i}"
done

# Load/activate env
source ~/bin/myconda
conda activate "${OPENFE_ENV}"

# Threading (still relevant for CPU-side tasks)
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# Optional: Some OpenMM installs respect these
export CUDA_VISIBLE_DEVICES="${SLURM_LOCALID:-0}"
                                           
mapfile -t JSONS < <(find "${JSON_DIR}" -maxdepth 1 -type f -name "*.json" | sort)

N_JSON="${#JSONS[@]}"
if [[ "${N_JSON}" -eq 0 ]]; then
  echo "No JSON files found in ${JSON_DIR}"
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID}"
JSON_INDEX=$(( TASK_ID / N_REP ))
REP_INDEX=$(( TASK_ID % N_REP ))

if [[ "${JSON_INDEX}" -ge "${N_JSON}" ]]; then
  echo "JSON_INDEX ${JSON_INDEX} out of range (N_JSON=${N_JSON}). Check --array range."
  exit 1
fi

JSON_FILE="${JSONS[JSON_INDEX]}"
BASE="$(basename "${JSON_FILE}" .json)"

WORKDIR="${OUT_ROOT}/replicate_${REP_INDEX}/${BASE}"
RESULT_JSON="${OUT_ROOT}/replicate_${REP_INDEX}/${BASE}.json"

mkdir -p "${WORKDIR}"

echo "Running:"
echo "  JSON_FILE   = ${JSON_FILE}"
echo "  WORKDIR     = ${WORKDIR}"
echo "  RESULT_JSON = ${RESULT_JSON}"
echo "  GPU         = 1"

openfe quickrun "${JSON_FILE}" -d "${WORKDIR}" -o "${RESULT_JSON}"                                                                                          
