//
// Created by lukas on 02.09.20.
//

#ifndef EXAMPLES_COMPLEX_POYLNOMIAL_MODEL_HPP
#define EXAMPLES_COMPLEX_POYLNOMIAL_MODEL_HPP

#include "../../include/simulation_header.hpp"
#include "execution/executer.hpp"


// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_complex_polynomial_model_complex_langevin()
{
    typedef std::complex<double> BasicType;
    typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParams;
    typedef ComplexLangevinUpdateParameters<ModelParams> MCMCUpdateParams;
    typedef SiteSimpleUpdateParameters UpdateDynamicsParams;
    typedef SiteParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "PolySiteModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(0.005);

    ModelParams model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "ComplexConfig"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(500, 1000000, 10000, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            site_parameters, execution_parameters, rel_config_path, rel_data_path);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

// ### Working example ###

void example_complex_polynomial_model_cobrid_monte_carlo()
{
    typedef double BasicType;
    typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParams;
    typedef CobridMonteCarloUpdateParameters<BasicType, ModelParams> MCMCUpdateParams;
    typedef GlobalSiteUpdateParameters UpdateDynamicsParams;
    typedef SiteParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "PolySiteModelCobridMonteCarlo";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(json {
            {"dt", 0.02},
            {"n", 100}
    });

    ModelParams model_parameters({1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(20, 100000, 10000, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            site_parameters, execution_parameters, rel_config_path, rel_data_path);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

// ### Working example ###

void example_complex_polynomial_model_cobrid_imag_monte_carlo()
{
    typedef double BasicType;
    typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParams;
    typedef CobridImagMonteCarloUpdateParameters<BasicType, ModelParams> MCMCUpdateParams;
    typedef GlobalSiteUpdateParameters UpdateDynamicsParams;
    typedef SiteParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "PolySiteModelCobridImagMonteCarlo";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(json {
            {"dt", 0.02},
            {"n", 100}
    });

    ModelParams model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "Config"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
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
