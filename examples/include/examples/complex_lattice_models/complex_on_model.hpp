//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONSEXAMPLES_COMPLEX_ON_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONSEXAMPLES_COMPLEX_ON_MODEL_HPP


#include "../../simulation_header.hpp"
#include "execution/executer.hpp"

void example_complex_on_model_complex_langevin()
{
    typedef lm_impl::link::ON<std::complex<double>, 4> BasicType;
    typedef lm_impl::lattice_system::ComplexONModelParameters ModelParams;
    typedef lm_impl::mcmc_update::ComplexLangevinUpdateONParameters<ModelParams> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::ParallelUpdateParameters UpdateDynamicsParams; // /*WithAdpativeStepsize*/
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "ComplexONModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.001);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0, 0.1, 0.1);

    UpdateDynamicsParams update_dynamics_parameters; // (2000);

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Mean", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 10, 0, {}, // optional additional measures
                                         {}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, rel_data_path, "model_params", "beta", 0.1, 0.6, 2);
    simparams.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}


#endif //LATTICEMODELIMPLEMENTATIONSEXAMPLES_COMPLEX_ON_MODEL_HPP
