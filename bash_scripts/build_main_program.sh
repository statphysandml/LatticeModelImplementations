#!/bin/bash

# Project name
read -p "Enter main project name: " project_name

if test -z $project_name; then
	echo "Project name cannot be empty."
	exit 1
fi

project_path="../"

if test -d "$project_path/include"; then
	echo "Main project already build. Only the cmake lists file is rebuild."

	mkdir -p "$project_path/cmake"
	# Generate CMakeLists.txt file
	source "${path_to_lattice_model_implementations}bash_scripts/write_cmake_lists_file_main_project.sh"
else
	# Generate folders for main project
	mkdir "$project_path/include"
	mkdir "$project_path/src"
	mkdir "$project_path/jupyter_notebooks"

	include_path="$project_path/include/"
	src_path="$project_path/src/"

	source "${path_to_lattice_model_implementations}bash_scripts/build_structure.sh"
fi
