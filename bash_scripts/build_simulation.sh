#!/bin/bash

# MCMC Simulation lib
if [ -z ${path_to_mcmc_simulation_lib+x} ]; then
  parent_dir="$(dirname -- "$(readlink -f -- "build_project.sh")")"
  path_to_lattice_model_implementations="$(dirname "$parent_dir")"
fi

path_to_base_lib=${path_to_lattice_model_implementations}

# Submodules
path_to_mcmc_simulation_lib="${path_to_lattice_model_implementations}/external_submodules/MCMCSimulationLib/"
path_to_param_helper="${path_to_mcmc_simulation_lib}/external_submodules/ParamHelper/"

# Determine project_name, project_path and project_type
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_simulation_builder.sh"

# Build the project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_main_builder.sh"

# Compile the sample project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_compiling.sh"