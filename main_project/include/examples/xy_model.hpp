//
// Created by lukas on 31.08.20.
//

#ifndef EXAMPLES_XY_MODEL_HPP
#define EXAMPLES_XY_MODEL_HPP

#include "../../include/simulation_header.hpp"
#include "execution/executer.hpp"

#include "../../../include/lattice/lattice_models/xy_model.hpp"

// ./LatticeModelImplementations mode_type config_file data_root_dir

// expectation_value XYModel

// ### Working example ###

void example_xy_model_metropolis()
{
    typedef double BasicType;
    typedef XYModelParameters<GaussianSampler> ModelParams;
    typedef MetropolisUpdateParameters<ModelParams> MCMCUpdateParams;
    typedef SequentialUpdateParameters UpdateDynamicsParams;
    typedef LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "XYModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters;

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"XYMagnetization"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 10000, 100, {}, // optional additional measures
                                         {}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            lattice_parameters, execution_parameters, rel_config_path, rel_data_path, "model_params", "beta", 0.05, 1.55, 10);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

// expectation_value XYModelHMC

// ### Working example ###

void example_xy_model_hmc_algorithm()
{
    typedef double BasicType;
    typedef XYModelParameters<GaussianSampler> ModelParams;
    typedef HybridMonteCarloUpdateParameters<BasicType , ModelParams > MCMCUpdateParams;
    typedef GlobalLatticeUpdateParameters UpdateDynamicsParams;
    typedef LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "XYModelHMC";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(json {
            {"dt", 0.02},
            {"n", 20}
    });

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"XYMagnetization"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(1, 10000, 100, {}, // optional additional measures
                                         {}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            lattice_parameters, execution_parameters, rel_config_path, rel_data_path, "model_params", "beta", 0.05, 1.55, 10);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}


#endif //EXAMPLES_XY_MODEL_HPP
