git submodule update --init --recursive

# Build ParamHelper
cd ../external_submodules/MCMCSimulationLib/external_submodules/ParamHelper/build
bash build.sh

# Build MCMCSimulationLib
cd ../../../build/
source build.sh

# Navigate back to build directory
cd ../../../build/