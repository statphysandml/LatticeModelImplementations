path_to_python3="~/.miniconda3/envs/flowequation/" # (optional)
virtual_env="flowequation" # (optional)
python_version="3.7" # (optional)
path_to_conda_activate="~/.miniconda3/bin/activate" # (optional)

# Only necessary, if GPU is used
# Laptop
nvcc_flag_gencode_arch=compute_75
nvcc_flag_gencode_code=sm_75

# GPU Cluster
# nvcc_flag_gencode_arch=compute_60
# nvcc_flag_gencode_code=sm_60
