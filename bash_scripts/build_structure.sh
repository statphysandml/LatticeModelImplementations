#!/bin/bash

# Generate folder structure
mkdir "$project_path/cmake"
mkdir "$project_path/debug"
mkdir "$project_path/release"

# Generate CMakeLists.txt file
source "${path_to_lattice_model_implementations}bash_scripts/write_cmake_lists_file_main_project.sh"

# Generate main.cpp file
source "${path_to_lattice_model_implementations}bash_scripts/write_main_cpp_file.sh"

# Generate simulation_header.hpp file
source "${path_to_lattice_model_implementations}bash_scripts/write_simulation_header_hpp_file.sh"

# Generate config_h_in file
source "${path_to_lattice_model_implementations}bash_scripts/write_config_h_in_file.sh"