#!/bin/bash

# Project name
read -p "Enter project name: " project_name

if test -z $project_name; then
	echo "Project name cannot be empty."
	exit 1
fi

project_path="../projects/${project_name}/"

if test -d "$project_path"; then
	echo "Project already exists."
  echo "Main project already build. Only the cmake lists file is rebuild."

	mkdir -p "$project_path/cmake"
	mkdir -p "$project_path/debug"
	mkdir -p "$project_path/release"
	# Generate CMakeLists.txt file
	source "${path_to_lattice_model_implementations}bash_scripts/write_cmake_lists_file_main_project.sh"
else

  mkdir -p "../projects/"
  mkdir $project_path

  include_path=$project_path
  src_path=$project_path

  source "${path_to_lattice_model_implementations}bash_scripts/build_structure.sh"

fi