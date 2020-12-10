//
// Created by lukas on 07.08.20.
//

#ifndef EXAMPLES_COMPLEX_XY_MODEL_HPP
#define EXAMPLES_COMPLEX_XY_MODEL_HPP


#include "../../simulation_header.hpp"
#include "execution/executer.hpp"

void example_complex_xy_model_complex_langevin()
{
    typedef double BasicType;
    typedef lm_impl::lattice_system::ComplexXYModelParameters ModelParams;
    typedef lm_impl::mcmc_update::ComplexLangevinUpdateParameters<ModelParams, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::ParallelUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "ComplexXYModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.001);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Mean"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 10000, 100, {}, // optional additional measures
                                         {}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, rel_data_path, "model_params", "beta", 0.05, 1.55, 10);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

#endif //EXAMPLES_COMPLEX_XY_MODEL_HPP
