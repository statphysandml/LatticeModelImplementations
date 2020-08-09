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

    Py_Initialize();
    PyRun_SimpleString("import sys\n" "import os");
    PyRun_SimpleString("sys.path.append( os.path.dirname(os.getcwd()) + '/')");
    PyRun_SimpleString("sys.path.append( os.path.dirname(os.getcwd()) + '/python_scripts/')");
    PyRun_SimpleString("sys.path.append( os.path.dirname(os.getcwd()) + '/python_scripts/plotting_environment/')");
    PyRun_SimpleString("print('Running python in ' + os.path.dirname(os.getcwd()) + '/python_scripts/')");

    if(argc > 1)
        run_from_file<SystemBaseParams>(argc, argv);
    else
        custom_main();

    Py_Finalize();
}

void custom_main()
{
    std::vector<int> dimensions {4, 4};

    // MetropolisUpdateParameters<IsingModelParameters> update_parameters;
    HybridMonteCarloUpdateParameters<IsingModelParameters> update_parameters(json {
            {"epsilon", 0.02},
            {"n", 100}
    });

    double beta = 0.4;
    IsingModelParameters ising_model_parameters(json {
            {"beta", beta},
            {"J", {1.0, 0.0}},
            {"h", {0.0, 0.0}},
            {"eps", 0.1}
    });

    // typedef LatticeParameters< int, IsingModelParameters, MetropolisUpdateParameters<IsingModelParameters>, SequentialUpdateParameters> LatticeParams;
    typedef LatticeParameters< int, IsingModelParameters, HybridMonteCarloUpdateParameters<IsingModelParameters>, GlobalLatticeUpdateParameters> LatticeParams;

    LatticeParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"lattice_update_type", "sequential"},
                    {"thermalization_steps", 20000},
                    {"measures", {"Mean", "AbsMean", "Std", "Energy", "Config"}},
                    {"model", ising_model_parameters.get_json()},
                    {"update", update_parameters.get_json()}},
                    "None"
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 1000, 1000, {}, // optional additional measures
                                         {"Mean", "Std", "Energy"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = SimulationParameters< LatticeParams, ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, "/data/IsingModel/", "model", "beta", 0.1, 0.625, 21);

    Simulation<LatticeParams, ExpectationValueParameters > sim(simparams);

    sim.run();
}

// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel
