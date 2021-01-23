path_to_python3="~/.miniconda3/envs/virtual_env/" # (optional)
virtual_env="virtual_env" # (optional)
python_version="3.7" # (optional)
path_to_conda_activate="~/.miniconda3/bin/activate" # (optional)

# Optional

# path_to_boost="/opt/boost_1_70_0" # can be defined if CMake cannot find boost

# Only necessary, if GPU is used - adapt this to your architecture
# path_to_cuda="/opt/cuda-10.1" # can be defined if CMake cannot find cuda
# nvcc_flag_gencode_arch=compute_75
# nvcc_flag_gencode_code=sm_75