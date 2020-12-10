//
// Created by lukas on 11.09.20.
//

#ifndef MAIN_CPP
#define MAIN_CPP

#include <Python.h>

#include "mcmc_simulation/header.hpp"

#include "../include/lattice_model_impl/update_dynamics/update_dynamics_header.hpp"
#include "../include/lattice_model_impl/mcmc_update/mcmc_update_header.hpp"

#include "../include/lattice_model_impl/site/site_header.hpp"
#include "../include/lattice_model_impl/lattice/lattice_header.hpp"
#include "../include/lattice_model_impl/link_lattice/link_lattice_header.hpp"

#include "execution/executer.hpp"

// #include "mcmc_simulation/simulation.hpp"

// #include "glog/logging.h"

#endif

// Todo next: Rename lattice_update_formalism to site_update_formalism and write new base class called lattice_update_formalism to be able to perform different kinds of update, these will also contain the necessary static asserts
// - Rewrite factory method

void custom_main();

/*

template<typename T>
using can_update = typename detail::is_updateable<T, MetropolisUpdateParameters<IsingModelParameters, mcmc::sampler::GaussianSampler>::MCMCUpdate>;


template <class T>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
foo(T t) {
    std::cout << "foo<arithmetic T>\n";
    return t;
}

struct meow
{
    meow() {}
}; **/

// ToDo: Generate own structs for datatype bases?

int main(int argc, char **argv) {
    // google::InitGoogleLogging(argv[0]);

    // std::cout << (can_update<double>::value) << std::endl;
    // std::cout << (can_update<int>::value) << std::endl;

    // std::cout << foo(5) << std::endl;
//     std::cout << foo(meow()) << std::endl;

#ifdef PYTHON
    mcmc::execution::initialize_python();
#endif

    custom_main();

#ifdef PYTHON
    mcmc::execution::finalize_python();
#endif
}

// - Introduce a default update_dynamics for sites - not for lattices!


void custom_main()
{}

// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel
