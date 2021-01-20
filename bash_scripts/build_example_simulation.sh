# Build simulation example (from existing main.cpp and simulation_header.hpp

project_name=$simulation_project_name
project_path="./simulations/${simulation_project_name}/"
mkdir -p $project_path
project_path="$(cd "$project_path" && pwd -P)"
project_type="simulation"

# Build the project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_main_builder.sh"

# Compile the sample project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_compiling.sh"

# Navigate back to examples/
cd ../../
