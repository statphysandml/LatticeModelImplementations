#ifndef MAIN_CPP
#define MAIN_CPP

#include "../include/config.h"
#include "../include/simulation_header.hpp"
#include "execution/executer.hpp"

#endif

void custom_main();

int main(int argc, char **argv) {
    param_helper::fs::prfs::set_relative_path_to_project_root_dir("/../");

#ifdef PYTHON
    mcmc::execution::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    if(argc > 1)
        mcmc::execution::run_from_file<from_file_simulation::SystemBaseParams>(argc, argv);
    else
        custom_main();

#ifdef PYTHON
    mcmc::execution::finalize_python();
#endif
    return 0;
}

void custom_main()
{}

// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel
