#!/bin/bash

# For an entire project

# MCMC Simulation lib
parent_dir="$(dirname -- "$(readlink -f -- "build_project.sh")")"
path_to_lattice_model_implementations="$(dirname "$parent_dir")"

path_to_base_lib=${path_to_lattice_model_implementations}

# Submodules
path_to_mcmc_simulation_lib="${path_to_lattice_model_implementations}/external_submodules/MCMCSimulationLib/"
path_to_param_helper="${path_to_mcmc_simulation_lib}/external_submodules/ParamHelper/"

# Determine project_name, project_path and project_name
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_project_builder.sh"

# Build the project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_main_builder.sh"

# Compile the sample project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_compiling.sh"

# Add bash_script for possibility to generate simulations
mkdir -p "$project_path/bash_scripts/"
source "${path_to_lattice_model_implementations}/bash_scripts/write_project_related_build_simulation_sh_file.sh"