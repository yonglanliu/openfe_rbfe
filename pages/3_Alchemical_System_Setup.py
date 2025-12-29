from __future__ import annotations

import streamlit as st
from pathlib import Path
from openfe_app.config import ProjectPaths
import pandas as pd
from openfe_app.alchemical_setup_ui import alchemical_setup_tab
import os


if "project_root" not in st.session_state:
    st.warning("No project selected yet. Go to the Project page and create/select one.")
    st.stop()
paths = ProjectPaths(Path(st.session_state["project_root"]))
st.session_state["project_root"] = str(paths.root)

st.sidebar.write("Folders:")
st.sidebar.code(
    f"inputs:   {paths.inputs}\n"
    f"prepared:  {paths.planned}\n"
    f"results:  {paths.results}\n"
    f"analysis: {paths.analysis}"
)

# ---------- Tabs ----------
tab_plan, tab_run= st.tabs(
    ["1) System Setup", "2) Configuration"]
)
with tab_plan:
    #-----------------------------------------
    sdf_file = paths.inputs / "ligands.sdf"
    pdb_file = paths.inputs / "protein.pdb"

    if not (sdf_file.exists() and pdb_file.exists()):
        st.warning("Upload ligands.sdf and protein.pdb in Inputs first.")
        st.stop()

    # Call exactly ONCE, with a unique prefix
    alchemical_setup_tab(
        inputs_sdf=sdf_file,
        protein_pdb=pdb_file,
        planned_root=paths.planned
    )

    st.divider()
with tab_run:
    st.markdown("### Use CIL to run the Alchemical FEP simulation")
    c1, c2 = st.columns([2, 2])
    with c1:
        n_replicates = int(st.text_input(
            "Number of replicates to run", 
            value=3
        ))
    st.info(f"""Important: If you choose 3 parallel schedule, choose 1 here.
            If you choose 2 or more here, you will run 3 * 2 = 6 or more simulations for a system.""")
    with c2:
        env_name = st.text_input(
            "Environment Name", 
            value="openfe"
        )
    st.divider()
    st.markdown("### Configuration for HPC  Slurm Batch Jobs ")
    c1, c2, c3, c4, c5, c6, c7 = st.columns([1,1,1,1,1,1,1])
    with c1:
        job_name = st.text_input(
            "Job Name",
            value="rbfe"
        )
    with c2:
        platform = st.selectbox(
            "Platform",
            ["gpu", "cpu"],
            index=0
        )

        if platform == "gpu":
            cc1, cc2 = st.columns([2,2])
            with cc1:
                n_gpu = st.number_input(
                    "Number of GPU",
                    min_value=1,
                    value=1,
                    step = 1
                )
            with cc2:
                gpu_type = st.text_input("GPU Type",
                                    value="v100x")
    with c3:
        mem = int(st.text_input("Memory (G)", value=16))
    with c4:
        cpus_per_task = int(st.text_input("cpus-per-task", value=16))
    with c5:
        wall_time = int(st.text_input(
            "Wall Time (h)",
            value=120
        ))
    with c6:
        n_json_files = st.number_input(
            "Number of JSON Files to run",
            min_value=1,
            value=18,
            step=1
        )
        total_tasks = int(n_json_files) * int(n_replicates)
    with c7:
        n_jobs = st.number_input("Total Jobs", value=total_tasks, disabled=True)
    st.divider()

    if platform == "gpu":
        st.markdown("#### GPU configuration")
        st.session_state["gpu config"] = f"""#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --partition={platform} 
#SBATCH --time={wall_time}:00:00
#SBATCH --gres=gpu:{gpu_type}:{n_gpu}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}G
#SBATCH --array=0-{n_jobs - 1}
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err     

# -------- User settings --------
OPENFE_ENV=\"{env_name}\"               # conda env name
JSON_DIR=\"{paths.planned}\"     # folder containing 18 *.json files
OUT_ROOT=\"{paths.results}\"        # where replicate_* folders will be created
N_REP={n_replicates}                                # replicates per transformation
# --------------------------------
"""
        st.session_state["gpu config"] += ("""
mkdir -p logs
for ((i=0; i<${N_REP}; i++)); do
    mkdir -p "${OUT_ROOT}/replicate_${i}"
done

# Load/activate env
source ~/.bashrc
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
""")

        script= st.session_state["gpu config"]
        st.code(script, language="bash")
        st.download_button(
        "Download Configuration",
        data = script,
        file_name="gpu_config.sh",
        mime="text/x-shellscript")
#-------------------------------
#      CPU configuration
#-------------------------------
    if platform == "cpu":
        st.markdown("#### CPU configuration")
        st.session_state["cpu config"] = f"""#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --partition={platform}                # <-- change to your CPU partition
#SBATCH --time={wall_time}:00:00
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}G
#SBATCH --array=0-{n_jobs - 1}
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
     

# -------- User settings --------
OPENFE_ENV=\"{env_name}\"               # conda env name
JSON_DIR=\"{paths.planned}\"     # folder containing 18 *.json files
OUT_ROOT=\"{paths.results}\"        # where replicate_* folders will be created
N_REP={n_replicates}                                # replicates per transformation
# --------------------------------
"""
        st.session_state["cpu config"] += ("""                                          
mkdir -p logs
                                           
for ((i=0; i<${N_REP}; i++)); do
    mkdir -p "${OUT_ROOT}/replicate_${i}"
done

# Load/activate env (edit to your cluster style)
# module load anaconda/...
source ~/.bashrc
conda activate "${OPENFE_ENV}"
                                           
# CPU threading knobs (OpenMM uses these)
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK}" 

# Collect JSON files (sorted, stable ordering)
mapfile -t JSONS < <(find "${JSON_DIR}" -maxdepth 1 -type f -name "*.json" | sort) 
N_JSON="${#JSONS[@]}"
if [[ "${N_JSON}" -eq 0 ]]; then
  echo "No JSON files found in ${JSON_DIR}"
  exit 1
fi

# Decode array index -> (json_index, replicate_index)
TASK_ID="${SLURM_ARRAY_TASK_ID}"
JSON_INDEX=$(( TASK_ID / N_REP ))
REP_INDEX=$(( TASK_ID % N_REP ))       
                                           
# Decode array index -> (json_index, replicate_index)
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
echo "  CPUS        = ${SLURM_CPUS_PER_TASK}"

# Run the transformation
openfe quickrun "${JSON_FILE}" -d "${WORKDIR}" -o "${RESULT_JSON}"                                                                                            
""")

        script= st.session_state["cpu config"]
        st.code(script, language="bash")
        st.download_button(
        "Download Configuration",
        data = script,
        file_name="cpu_config.sh",
        mime="text/x-shellscript")
    st.divider()