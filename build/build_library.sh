#!/bin/bash

# Submodules
path_to_mcmc_simulation_lib="./external_submodules/MCMCSimulationLib/"
path_to_param_helper="${path_to_mcmc_simulation_lib}external_submodules/ParamHelper/"

# Remove lib directory if already existing to avoid errors when changing parameters
if [ -d "../lib" ]; then
  echo "Remove lib folder to avoid errors due to changing parameters."
  rm -r "../lib"
fi

# To extract necessary parameters for the cmake file
source "${path_to_config}/config.sh"

# Generate generate_cmake_lists_txt_file
source generate_cmake_lists_txt_file.sh

cd ../
mkdir -p lib
cd lib

cmake ..
make -j9

cd ../build/
