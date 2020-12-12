#ifndef MAIN_CPP
#define MAIN_CPP

#include "../include/simulation_header.hpp"
#include "execution/executer.hpp"

#endif

void custom_main();

int main(int argc, char **argv) {
    param_helper::fs::prfs::set_relative_path_to_project_root_dir("/../");

#ifdef PYTHON
    mcmc::execution::initialize_python();
#endif

    if(argc > 1)
        mcmc::execution::run_from_file<from_file_simulation::SystemBaseParams>(argc, argv);
    else
        custom_main();

#ifdef PYTHON
    mcmc::execution::finalize_python();
#endif
}

/* #include "../include/examples/complex_lattice_models/complex_anharmonic_oscillator.hpp"
#include "../include/examples/lattice_models/ising_model.hpp"
#include "../include/examples/lattice_models/xy_model.hpp"
#include "../include/examples/complex_site_models/complex_cubic_model.hpp"
#include "../include/examples/complex_site_models/complex_polynomial_model.hpp"
#include "../include/examples/complex_site_models/complex_scalar_gaussian_model.hpp"
#include "../include/examples/complex_lattice_models/complex_xy_model.hpp" */

// #include "../include/examples/lattice_models/ising_model_metropolis.hpp"
// #include "../include/examples/lattice_models/xy_model_hybrid_monte_carlo.hpp"
#include "../include/examples/lattice_models/u_one_model.hpp"
#include "../include/examples/lattice_models/su_two_model.hpp"
#include "../include/examples/lattice_models/on_model_metropolis.hpp"

void custom_main()
{
    // example_ising_model_metropolis();
    // example_xy_model_hmc_algorithm();
    example_u_one_model_metropolis();
    example_su_two_model_metropolis();
    example_on_model_metropolis();

    // example_ising_full_simulation();
    // example_xy_model_metropolis();
    // example_xy_model_hmc_algorithm();

    // example_complex_cubic_model_complex_langevin();
    // example_complex_scalar_gaussian_model();
    // example_complex_polynomial_model_complex_langevin();
    // example_complex_polynomial_model_cobrid_monte_carlo();

    // example_complex_polynomial_model_complex_langevin_correlation_time();

    /* example_complex_xy_model_complex_langevin(); ToDo: Check correct action and drift term definitions, implement adaptive stepsize!
    example_complex_xy_model_config_comparison_complex_langevin(); */

    // example_complex_anharmonic_oscillator_cobrid_monte_carlo();
}

// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel
