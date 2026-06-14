#!/usr/bin/env bash

source ~/bin/myconda
conda activate openfe_m

module purge
module load CUDA/11.8.0
module load cuDNN/8.9.2/CUDA-11

CONF_DIR="simulation/septop/conf"
CONF_YAML="${CONF_DIR}/septop_conf.yaml"

WORKDIR=$(yq -r '.workdir' "${CONF_YAML}")
RESDIR=$(yq -r '.result_dir' "${CONF_YAML}")
JSONDIR=$(yq -r '.json_dir' "${CONF_YAML}")

mkdir -p "${WORKDIR}" "${RESDIR}" "${JSONDIR}"

CMD=(
    python script/run_septop.py
    --conf_file "${CONF_YAML}"
)

printf '%q ' "${CMD[@]}"
echo
"${CMD[@]}"
