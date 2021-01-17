//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONSEXAMPLES_COMPLEX_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONSEXAMPLES_COMPLEX_POLYNOMIAL_MODEL_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"

// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_complex_polynomial_model_complex_langevin()
{
    typedef std::complex<double> BasicType;
    typedef lm_impl::site_system::ComplexPolynomialModelParameters ModelParams;
    typedef lm_impl::mcmc_update::ComplexLangevinUpdateParameters<ModelParams, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::UpdateWithAdpativeStepsizeParameters UpdateDynamicsParams;
    typedef lm_impl::site_system::SiteParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "ComplexPolynomialModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(0.01);

    ModelParams model_parameters(json {
        {"lambda_real", 1.0},
        {"lambda_imag", 0.0},
        {"sigma_real", 1.0},
        {"sigma_imag", 1.0},
        {"h_real", 1.0},
        {"h_imag", 0.0}}
    )

    UpdateDynamicsParams update_dynamics_parameters(2000);

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "ComplexConfig"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 100000, 0, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            site_parameters, execution_parameters, rel_data_path);
    simparams.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

#endif //LATTICEMODELIMPLEMENTATIONSEXAMPLES_COMPLEX_POLYNOMIAL_MODEL_HPP
