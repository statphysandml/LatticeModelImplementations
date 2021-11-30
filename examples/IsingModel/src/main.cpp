#ifndef MAIN_CPP
#define MAIN_CPP

#include "../include/IsingModel/config.h"
#include "../include/IsingModel/simulation_header.hpp"

#endif

void custom_main();

int main(int argc, char **argv) {
    param_helper::fs::prfs::set_relative_path_to_project_root_dir("../");

    // Initialization - Only needed for GPU and CPU runs
    mcmc::execution::initialize_executer_params(PROJECT_NAME, CLUSTER_MODE);

#ifdef RUN_WITH_PYTHON_BACKEND
    mcmc::execution::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    // The function of the first if-condition is only called when an actual simulation takes place based on arguments

    if(argc > 1)
    {
        mcmc::execution::run_from_file<typename from_file_simulation::SystemBaseParams>(argc, argv);
    }
    else
        // Helpful for a preparation of the simulation or immediate execution (on cpu/gpu/locally, testing/running directly)
        custom_main();

#ifdef RUN_WITH_PYTHON_BACKEND
    mcmc::execution::finalize_python();
#endif
    return 0;
}

void custom_main()
{
    typedef lm_impl::lattice_system::IsingModelParameters ModelParams;
    
    typedef double BasicType;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::lattice_system::IsingModelSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "IsingModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";
    std::string correlation_time_results_path = "/results/" + model_name + "/";

    std::vector<int> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.1);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Config", "Mean", "SecondMoment"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    /* mcmc::execution::EquilibriateParameters equilibriate_parameters(5000, 100, {"Mean"}); // Meausures which will be evaluated in terms of mean and error evaluation
    // equilibriate_parameters.write_to_file(rel_data_path);

    auto simulation_params_equilibriate = mcmc::simulation::SimulationParameters< SystemBaseParams, mcmc::execution::EquilibriateParameters >::generate_simulation(
            lattice_parameters, equilibriate_parameters, rel_data_path, "model_params", "beta", 0.1, 0.6, 5);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::EquilibriateParameters::name(), model_name); */

    mcmc::execution::CorrelationTimeParameters correlation_time_parameters(5000, 100, 10000, {"Mean"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simulation_params_correlation_time = mcmc::simulation::SimulationParameters< SystemBaseParams , mcmc::execution::CorrelationTimeParameters >::generate_simulation(
            lattice_parameters, correlation_time_parameters, rel_data_path, "model_params", "beta", 0.1, 0.7, 25);
    simulation_params_correlation_time.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::CorrelationTimeParameters::name(), model_name);

    mcmc::execution::ExpectationValueParameters expectation_value_parameters(correlation_time_results_path, 10000, 10000, {}, // optional additional measures
                                                                             {"AbsMean", "Energy"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simulation_params_expectation_value = mcmc::simulation::SimulationParameters< SystemBaseParams , mcmc::execution::ExpectationValueParameters >::generate_simulation(
            lattice_parameters, expectation_value_parameters, rel_data_path, "model_params", "beta", 0.1, 0.7, 25);
    simulation_params_expectation_value.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::ExpectationValueParameters::name(), model_name);
}

// Rerun the simulation with ./IsingModel expectation_value IsingModelMetropolis
