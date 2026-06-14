#!/bin/bash

#SBATCH --job-name=abfe
#SBATCH --partition=gpu 
#SBATCH --time=240:00:00
#SBATCH --gres=gpu:v100x:1,lscratch:100
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --array=0-2


# Load/activate env
source ~/bin/myconda
conda activate openfe_m

module purge
module load CUDA/11.8.0
module load cuDNN/8.9.2/CUDA-11


# -------- User settings --------
CONF_DIR="simulation/abfe/PDE4/conf"     
CONF_YAML="${CONF_DIR}/abfe_conf.yaml"

JSON_DIR=$(yq -r '.json_dir' "${CONF_YAML}")    
OUT_ROOT=$(yq -r '.result_dir' "${CONF_YAML}")  

echo "OUT_ROOT=${OUT_ROOT}"

N_REP=3                                # replicates per transformation


for ((i=0; i<${N_REP}; i++)); do
    mkdir -p "${OUT_ROOT}/replicate_${i}"
done


# Threading (still relevant for CPU-side tasks)
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# Optional: Some OpenMM installs respect these
export CUDA_VISIBLE_DEVICES="${SLURM_LOCALID:-0}"


# -------- Define if run json files in a JSON_DIR or run specific input_json_files -----
INPUT_JSON=$(yq -o=json '.input_json_files' "${CONF_YAML}")

# Case 1: not provided
INPUT_JSON=$(yq -r '.input_json_files // empty' "${CONF_YAML}")

if [[ -z "$INPUT_JSON" ]]; then
    echo "No json file specified. Will run all JSON files in ${JSON_DIR}"
    mapfile -t JSONS < <(find "${JSON_DIR}" -maxdepth 1 -type f -name "*.json" | sort)

elif [[ "$INPUT_JSON" == *.txt ]]; then
    echo "Reading JSON file list from $INPUT_JSON"
    if [[ ! -f "$INPUT_JSON" ]]; then
        echo "Error: list file not found: $INPUT_JSON"
        exit 1
    fi
    mapfile -t JSONS < "$INPUT_JSON"

elif [[ -f "$INPUT_JSON" ]]; then
    echo "Reading single JSON file: $INPUT_JSON"
    JSONS=("$INPUT_JSON")

else
    echo "Error: input_json_files is set but not a valid .txt or .json file: $INPUT_JSON"
    exit 1
fi


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

JSON_FILE="${JSONS[$JSON_INDEX]}"
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
