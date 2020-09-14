//
// Created by lukas on 10.09.20.
//

#ifndef EXAMPLES_COMPLEX_ANHARMONIC_OSCILLATOR_HPP
#define EXAMPLES_COMPLEX_ANHARMONIC_OSCILLATOR_HPP

#include "../../simulation_header.hpp"
#include "execution/executer.hpp"


void example_complex_anharmonic_oscillator_cobrid_monte_carlo()
{
    typedef double BasicType;
    typedef ComplexAnharmonicOscillatorParameters<GaussianSampler> ModelParams;
    typedef CobridMonteCarloUpdateParameters<BasicType, ModelParams> MCMCUpdateParams;
    typedef GlobalLatticeUpdateParameters UpdateDynamicsParameters;
    typedef LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParameters> SystemBaseParams;

    std::string model_name = "ComplexAnharmonicOscillator";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<int> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(json {
            {"dt", 0.002},
            {"n", 5}
    });

    ModelParams model_parameters(1.0, {1.0, 0.0}, {1.0, 1.0}, {6.0, 0.0}, 0.1);

    SystemBaseParams site_parameters(
            json {
                    {"dimensions", dimensions},
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

#endif //EXAMPLES_COMPLEX_ANHARMONIC_OSCILLATOR_HPP
