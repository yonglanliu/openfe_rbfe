# Load/activate env
source ~/bin/myconda
conda activate openfe_m

module purge
module load CUDA/11.8.0
module load cuDNN/8.9.2/CUDA-11

sdf_path="/simulation/rbfe/example/inputs"
input_sdf=${sdf_path}/ligands.sdf
output_sdf=${sdf_path}/charged_ligands.sdf

openfe charge-molecules -M ${input_sdf} -o ${output_sdf} -n 4 -s bash/charge_settings.yaml --overwrite-charges
