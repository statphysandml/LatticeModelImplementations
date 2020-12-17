//
// Created by lukas on 07.08.20.
//

#ifndef EXAMPLES_ISING_MODEL_METROPOLIS_HPP
#define EXAMPLES_ISING_MODEL_METROPOLIS_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"


// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel

// ### Working example ###

void example_ising_model_metropolis()
{
    typedef double BasicType;
    typedef lm_impl::lattice_system::IsingModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::lattice_system::IsingModelSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "IsingModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<int> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.1);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Config", "Mean", "Std"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 10000, 10000, {}, // optional additional measures
                                         {"AbsMean", "Energy"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, rel_data_path, "model_params", "beta", 0.1, 0.625, 21);
    simparams.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}


void example_ising_full_simulation()
{
    typedef double BasicType;
    typedef lm_impl::lattice_system::IsingModelParameters ModelParams;
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
                    {"measures", {}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    mcmc::execution::EquilibriateParameters equilibriate_parameters(5000, 100, {"Mean"}); // Meausures which will be evaluated in terms of mean and error evaluation
    // equilibriate_parameters.write_to_file(rel_data_path);

    auto simparams_equilibriate = mcmc::simulation::SimulationParameters< SystemBaseParams, mcmc::execution::EquilibriateParameters >::generate_simulation(
            lattice_parameters, equilibriate_parameters, rel_data_path, "model_params", "beta", 0.1, 0.6, 5);
    simparams_equilibriate.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::EquilibriateParameters::name(), model_name);

    mcmc::execution::CorrelationTimeParameters correlation_time_parameters(5000, 100, 10000, {"Mean"}); // Meausures which will be evaluated in terms of mean and error evaluation
    // correlation_time_parameters.write_to_file(rel_config_path);

    auto simparams_correlation_time = mcmc::simulation::SimulationParameters< SystemBaseParams , mcmc::execution::CorrelationTimeParameters >::generate_simulation(
            lattice_parameters, correlation_time_parameters, rel_data_path, "model_params", "beta", 0.1, 0.6, 5);
    simparams_correlation_time.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::CorrelationTimeParameters::name(), model_name);

    mcmc::execution::ExpectationValueParameters expectation_value_parameters(correlation_time_results_path, 10000, 10000, {"Mean", "Config", "Energy"}, // optional additional measures
                                         {"AbsMean"}); // Meausures which will be evaluated in terms of mean and error evaluation
    // expectation_value_parameters.write_to_file(rel_config_path);

    auto simparams_expectation_value = mcmc::simulation::SimulationParameters< SystemBaseParams , mcmc::execution::ExpectationValueParameters >::generate_simulation(
            lattice_parameters, expectation_value_parameters, rel_data_path, "model_params", "beta", 0.1, 0.6, 5);
    simparams_expectation_value.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::ExpectationValueParameters::name(), model_name);
}


#endif //EXAMPLES_ISING_MODEL_METROPOLIS_HPP
