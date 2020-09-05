#ifndef MAIN_CPP
#define MAIN_CPP


#include "../include/simulation_header.hpp"
#include "execution/executer.hpp"

// #include "mcmc_simulation/simulation.hpp"

#include "glog/logging.h"

#endif

// Todo next: Rename lattice_update_formalism to site_update_formalism and write new base class called lattice_update_formalism to be able to perform different kinds of update, these will also contain the necessary static asserts
// - Rewrite factory method

void custom_main();

template<typename T>
using can_update = typename detail::is_updateable<T, MetropolisUpdateParameters<IsingModelParameters>::UpdateFormalism>;


template <class T>
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
foo(T t) {
    std::cout << "foo<arithmetic T>\n";
    return t;
}

struct meow
{
    meow() {}
};

// ToDo: Generate own structs for datatype bases?

int main(int argc, char **argv) {
    google::InitGoogleLogging(argv[0]);

    std::cout << (can_update<double>::value) << std::endl;
    std::cout << (can_update<int>::value) << std::endl;

    std::cout << foo(5) << std::endl;
//     std::cout << foo(meow()) << std::endl;

#ifdef GPU
    /* readdy::util::thrust_integration::IntegrationWeights<double> integration_weights;
    integration_weights.set_weights();
    std::cout << "go" << std::endl;

    funct func;
    std::pair<double, double> result;
    for(auto i = 0; i < 100000; i++) {
        result = readdy::util::integration::integrate([](double x) { return std::pow(x, 2) * std::sin(x) * std::exp(-1.0) * std::sin(x) * std::exp(-1.0) * std::sin(x) * std::exp(-1.0); }, -1.3, 3.3);
    }
    std::cout << result.first << " " << result.second << std::endl;

    thrust_funct thrust_func_finite;
    std::pair<double, double> thrustresult_finite;
    readdy::util::thrust_finite_integration::integrate_wrapper<thrust_funct, double> integration_wrapper(thrust_func_finite, 100);
    for(auto i = 0; i < 100000; i++) {
        thrustresult_finite = integration_wrapper.integrate(-1.3, 3.3);
    }
    std::cout << thrustresult_finite.first << " " << thrustresult_finite.second << std::endl;

    thrust_funct thrust_func;
    std::pair<double, double> thrustresult;
    for(auto i = 0; i < 100000; i++) {
        readdy::util::thrust_integration::integrate_wrapper<thrust_funct, double> integration_wrapper(thrust_func, integration_weights, -1.3, 3.3);
        thrustresult = integration_wrapper.integrate();
        // thrustresult = readdy::util::thrust_integration::integrate(thrust_func, integration_weights, -1.3, 3.3);
    }
    std::cout << thrustresult.first << " " << thrustresult.second << std::endl; */
#endif
    // funct func;
    // auto result = readdy::util::integration::integrate([] (double x) { return std::pow(x, 2); }, -1.3, 3.3);
    // auto result = readdy::util::integration::integrate(func, -1.3, 3.3);

    initialize_python();

    if(argc > 1)
        run_from_file<from_file_simulation::SystemBaseParams>(argc, argv);
    else
        custom_main();

    finalize_python();
}

#include "../include/examples/ising_model.hpp"
#include "../include/examples/xy_model.hpp"
#include "../include/examples/complex_cubic_model.hpp"
#include "../include/examples/complex_polynomial_model.hpp"
#include "../include/examples/complex_scalar_gaussian_model.hpp"

void custom_main()
{
    // example_ising_model();
    // example_xy_model_metropolis();
    // example_xy_model_hmc_algorithm();

    // example_complex_cubic_model_complex_langevin();
    // example_complex_scalar_gaussian_model();
    example_complex_polynomial_model_complex_langevin();
    // example_complex_polynomial_model_cobrid_monte_carlo();
    example_complex_polynomial_model_cobrid_imag_monte_carlo();
}

// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel
