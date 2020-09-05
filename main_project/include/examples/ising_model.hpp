//
// Created by lukas on 07.08.20.
//

#ifndef EXAMPLES_ISING_MODEL_HPP
#define EXAMPLES_ISING_MODEL_HPP

#include "../../include/simulation_header.hpp"
#include "execution/executer.hpp"


// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value IsingModel

// ToDo: Introduce default sampler if sampling is not used!! Or remove completely, if not needed

// ### Working example ###

void example_ising_model_metropolis()
{
    typedef double BasicType;
    typedef IsingModelParameters ModelParams;
    typedef MetropolisUpdateParameters<ModelParams> MCMCUpdateParams;
    typedef SequentialUpdateParameters UpdateDynamicsParams;
    typedef LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "IsingModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<int> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters;

    double beta = 0.4;
    ModelParams model_parameters(beta, {1.0, 0.0}, {0.0, 0.0});

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Mean", "Std", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 10000, 10000, {}, // optional additional measures
                                         {"AbsMean", "Energy"}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            lattice_parameters, execution_parameters, rel_config_path, rel_data_path, "model_params", "beta", 0.1, 0.625, 21);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}


#endif //EXAMPLES_ISING_MODEL_HPP
