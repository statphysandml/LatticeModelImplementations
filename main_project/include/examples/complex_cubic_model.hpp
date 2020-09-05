//
// Created by lukas on 07.08.20.
//

#ifndef EXAMPLES_CUBIC_GAUSSIAN_MODEL_HPP
#define EXAMPLES_CUBIC_GAUSSIAN_MODEL_HPP

#include "../../include/simulation_header.hpp"
#include "execution/executer.hpp"

// ./LatticeModelImplementations mode_type config_file data_root_dir

// ### Working example ###

void example_complex_cubic_model_complex_langevin()
{
    typedef std::complex<double> BasicType;
    typedef ComplexCubicModelParameters<GaussianSampler> ModelParams;
    typedef ComplexLangevinUpdateParameters<ModelParams> MCMCUpdateParams;
    typedef SiteSimpleUpdateParameters UpdateDynamicsParams;
    typedef SiteParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "CubicSiteModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    MCMCUpdateParams mcmc_update_parameters(0.02);

    ModelParams model_parameters(0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams site_parameters(
            json {
                    {"measures", {"Mean", "ComplexConfig"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            "None"
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 10000, 100, {}, // optional additional measures
                                         {}); // Meausures which will be evaluated in terms of mean and error evaluation
    execution_parameters.write_to_file(rel_data_path);

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            site_parameters, execution_parameters, rel_config_path, rel_data_path);

    execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

#endif //EXAMPLES_CUBIC_GAUSSIAN_MODEL_HPP
