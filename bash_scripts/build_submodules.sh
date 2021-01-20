git submodule update --init --recursive

# Build ParamHelper
cd ../external_submodules/LatticeModelSimulationLib/build
source build.sh

# Navigate back to bash_scripts/
cd ../../../bash_scripts

echo "Path ${PWD}"