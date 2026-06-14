#!/usr/bin/env bash

# Load/activate env
source ~/bin/myconda
conda activate openfe_m

module purge
module load CUDA/11.8.0
module load cuDNN/8.9.2/CUDA-11

CONF_DIR="simulation/abfe/conf"
CONF_YAML="${CONF_DIR}/abfe_conf.yaml"

WORKDIR=$(yq -r '.workdir' "${CONF_YAML}")
RESDIR=$(yq -r '.result_dir' "${CONF_YAML}")
JSONDIR=$(yq -r '.json_dir' "${CONF_YAML}")

mkdir -p "${WORKDIR}" "${RESDIR}" "${JSONDIR}"

# Build command as array (safe for spaces)
CMD=(
    python script/run_abfe.py
    --conf_file "${CONF_YAML}"
)

# Run
echo "${CMD[@]}"
"${CMD[@]}"
