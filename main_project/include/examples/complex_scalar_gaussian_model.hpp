//
// Created by lukas on 04.09.20.
//

#ifndef EXAMPLES_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP
#define EXAMPLES_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP

#include "../../include/simulation_header.hpp"
#include "execution/executer.hpp"

// ### Working example ###

void example_complex_scalar_gaussian_model()
{
    typedef std::complex<double> BasicType;
    typedef ComplexScalarGaussianModelParameters<GaussianSampler> ModelParams;
    typedef ComplexLangevinUpdateParameters<ModelParams > MCMCUpdateParams;
    typedef SiteSimpleUpdateParameters UpdateDynamicsParams;
    typedef SiteParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "ComplexScalarGaussianModel";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(json {
            {"epsilon", 0.001}
    });

    // ModelParams model_parameters({0.5, -1.0}, 0.0, 0.1);
    ModelParams model_parameters({1.0, 1.0}, 0.0, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"measures", {"Mean", "ComplexConfig"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 100000, 10000, {}, // optional additional measures
                                         {"2ndMoment"}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            lattice_parameters, execution_parameters, rel_config_path, rel_data_path);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}


#endif //EXAMPLES_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP
