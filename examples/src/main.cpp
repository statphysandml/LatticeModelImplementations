#ifndef MAIN_CPP
#define MAIN_CPP

#include "../include/simulation_header.hpp"
#include "execution/executer.hpp"

#endif


void custom_main();

int main(int argc, char **argv) {

    initialize_python();

    if(argc > 1)
        run_from_file<from_file_simulation::SystemBaseParams>(argc, argv);
    else
        custom_main();

    finalize_python();
}

/* #include "../include/examples/complex_lattice_models/complex_anharmonic_oscillator.hpp"
#include "../include/examples/lattice_models/ising_model.hpp"
#include "../include/examples/lattice_models/xy_model.hpp"
#include "../include/examples/complex_site_models/complex_cubic_model.hpp"
#include "../include/examples/complex_site_models/complex_polynomial_model.hpp"
#include "../include/examples/complex_site_models/complex_scalar_gaussian_model.hpp"
#include "../include/examples/complex_lattice_models/complex_xy_model.hpp" */
#include "../include/examples/integration/integration.hpp"
#include "../include/examples/implicit_solver/implicit_solver.hpp"

// - Introduce a default update_dynamics for sites - not for lattices!


void custom_main()
{
    // integrate();
    examples::implicit_solver::run_implicit_integral_solvers();
#ifdef THRUST
    examples::implicit_solver::thrust_run_implicit_integral_solvers();
#endif

    // example_ising_model();
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
