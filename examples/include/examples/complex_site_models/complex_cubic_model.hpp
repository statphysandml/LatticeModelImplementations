//
// Created by lukas on 07.08.20.
//

#ifndef EXAMPLES_CUBIC_GAUSSIAN_MODEL_HPP
#define EXAMPLES_CUBIC_GAUSSIAN_MODEL_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"

// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_complex_cubic_model_complex_langevin()
{
    typedef std::complex<double> BasicType;
    typedef lm_impl::site_system::ComplexCubicModelParameters ModelParams;
    typedef lm_impl::mcmc_update::ComplexLangevinUpdateParameters<ModelParams, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
    typedef lm_impl::site_system::SiteParameters< BasicType, ModelParams, MCMCUpdateParams> SystemBaseParams;

    std::string model_name = "CubicSiteModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(0.01);

    ModelParams model_parameters;

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 100000, 1000, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            site_parameters, execution_parameters, rel_data_path);
    simparams.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

#endif //EXAMPLES_CUBIC_GAUSSIAN_MODEL_HPP
