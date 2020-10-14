path_to_python3="/remote/lin34/kades/.conda/envs/pytorchlocal3/"
path_to_param_helper="/remote/lin34/kades/ParamHelper/"
path_to_mcmc_simulation_lib="/remote/lin34/kades/MCMCSimulationLib/"

python_version="3.8"


# Optional
path_to_boost="/opt/boost_1_70_0"

# Only necessary, if GPU is used
# Laptop
# nvcc_flag_gencode_arch=compute_75
# nvcc_flag_gencode_code=sm_75

# GPU Cluster
nvcc_flag_gencode_arch=compute_60
nvcc_flag_gencode_code=sm_60
