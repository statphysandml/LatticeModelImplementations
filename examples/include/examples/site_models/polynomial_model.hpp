//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONSEXAMPLES_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONSEXAMPLES_POLYNOMIAL_MODEL_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"

// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_polynomial_model_langevin()
{
    typedef double BasicType;
    typedef lm_impl::site_system::PolynomialModelParameters ModelParams;
    typedef lm_impl::mcmc_update::LangevinUpdateParameters<ModelParams, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
    typedef lm_impl::site_system::SiteParameters< BasicType, ModelParams, MCMCUpdateParams> SystemBaseParams;

    std::string model_name = "PolynomialModelLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(0.01);

    ModelParams model_parameters(1.0, 1.0, 0.0);

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 100000, 0, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            site_parameters, execution_parameters, rel_data_path);
    simparams.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

#endif //LATTICEMODELIMPLEMENTATIONSEXAMPLES_POLYNOMIAL_MODEL_HPP
