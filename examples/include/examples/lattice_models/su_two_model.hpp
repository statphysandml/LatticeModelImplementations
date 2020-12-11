//
// Created by lukas on 07.08.20.
//

#ifndef EXAMPLES_SU_TWO_MODEL_HPP
#define EXAMPLES_SU_TWO_MODEL_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"


// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_su_two_model_metropolis()
{
    typedef SU2<double> BasicType;
    typedef lm_impl::link_lattice_system::SUTwoModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::link_lattice_system::SUTwoModelSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "SUTwoModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<int> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.01);

    double beta = 0.4;
    ModelParams model_parameters(beta);

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

#endif //EXAMPLES_SU_TWO_MODEL_HPP
