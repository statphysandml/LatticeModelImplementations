//
// Created by lukas on 02.09.20.
//

#ifndef EXAMPLES_COMPLEX_POYLNOMIAL_MODEL_HPP
#define EXAMPLES_COMPLEX_POYLNOMIAL_MODEL_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"


// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_complex_polynomial_model_cobrid_monte_carlo()
{
    typedef double BasicType;
    typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParams;
    typedef SingleSiteCobridMonteCarloUpdateParameters<BasicType, ModelParams, complex_gaussian_distribution_from_file> MCMCUpdateParams;
    typedef SiteParameters< BasicType, ModelParams, MCMCUpdateParams> SystemBaseParams;

    std::string model_name = "PolySiteModelCobridMonteCarlo";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(json {
            {"dt", 0.02},
            {"n", 100}
    });

    ModelParams model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}, 0.1);

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(20, 200000, 10000, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            site_parameters, execution_parameters, rel_config_path, rel_data_path);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

#endif //EXAMPLES_COMPLEX_POYLNOMIAL_MODEL_HPP
