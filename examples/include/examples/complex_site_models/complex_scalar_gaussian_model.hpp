//
// Created by lukas on 04.09.20.
//

#ifndef EXAMPLES_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP
#define EXAMPLES_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"

// ### Working example ###

void example_complex_scalar_gaussian_model()
{
    typedef std::complex<double> BasicType;
    typedef lm_impl::site_system::ComplexScalarGaussianModelParameters ModelParams;
    typedef lm_impl::mcmc_update::ComplexLangevinUpdateParameters<ModelParams, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
    typedef lm_impl::site_system::SiteParameters< BasicType, ModelParams, MCMCUpdateParams> SystemBaseParams;

    std::string model_name = "ComplexScalarGaussianModel";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(json {
            {"epsilon", 0.001}
    });

    ModelParams model_parameters({0.2, 0.2}, 0.0, 0.1);

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "ComplexConfig"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 100000, 10000, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            site_parameters, execution_parameters, rel_data_path);
    simparams.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);

    mcmc::execution::PlotSiteDistributionParameters site_distribution_parameters(
            json {{"xkey", "StateReal"},
                  {"ykey", "StateImag"},
                  {"rmin_x", -5.0},
                  {"rmax_x", 5.0},
                  {"rmin_y", -5.0},
                  {"rmax_y", 5.0}});
    site_distribution_parameters.write_to_file(rel_data_path);

    auto simparams_plot_site_distribution = mcmc::simulation::SimulationParameters<
            SystemBaseParams, mcmc::execution::PlotSiteDistributionParameters >::generate_simulation(
                    site_parameters, site_distribution_parameters, rel_data_path);
    simparams_plot_site_distribution.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (mcmc::execution::PlotSiteDistributionParameters::name(), model_name);
}


#endif //EXAMPLES_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP
