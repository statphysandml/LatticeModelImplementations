#!/bin/bash

# Examples

# Set config.sh file directory
path_to_config="$(dirname -- "$(readlink -f -- "build_examples.sh")")/"
echo "Path to config ${path_to_config}, ${path_to_base_lib}"

# Build Submodules
source build_submodules.sh

# Load config files for generation of CMakeLists.txt files
# The first line is only necessary if submodules are not build a priori
source config.sh
source project_config.sh

# LatticeModelImplementations path
parent_dir="$(dirname -- "$(readlink -f -- "build_examples.sh")")"
path_to_lattice_model_implementations="$(dirname "$parent_dir")"

# Submodules
path_to_base_lib="${path_to_lattice_model_implementations}/external_submodules/LatticeModelSimulationLib/"
path_to_lattice_simulation_lib=${path_to_base_lib}
path_to_mcmc_simulation_lib="${path_to_base_lib}external_submodules/MCMCSimulationLib/"
path_to_param_helper="${path_to_mcmc_simulation_lib}external_submodules/ParamHelper/"

project_name="LatticeModelImplementations"
project_path="./../"
project_path="$(cd "$project_path" && pwd -P)"
project_type="project"



# Build the project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_main_builder.sh"

# Compile the sample project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_compiling.sh"

# Build example simulations

simulation_project_name="ComplexONModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="ComplexPolynomialModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="ComplexXYModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="IsingModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="ONModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="PolynomialModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="SU2Model"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="U1Model"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

simulation_project_name="XYModel"
source "${path_to_lattice_model_implementations}/bash_scripts/build_example_simulation.sh"

# simulation_project_name=""
# source build_example_simulation.sh
