#ifndef MAIN_CPP
#define MAIN_CPP

#include "../include/ComplexXYModel/config.h"
#include "../include/ComplexXYModel/simulation_header.hpp"

#endif

void custom_main();

int main(int argc, char **argv) {
    param_helper::fs::prfs::set_relative_path_to_project_root_dir("../");

    // Initialization - Only needed for GPU and CPU runs
    mcmc::execution::initialize_executer_params(PROJECT_NAME, CLUSTER_MODE);

#ifdef RUN_WITH_PYTHON_BACKEND
    mcmc::execution::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    // The function of the first if-condition is only called when an actual simulation takes place based on arguments

    if(argc > 1)
    {
        mcmc::execution::run_from_file<typename from_file_simulation::SystemBaseParams>(argc, argv);
    }
    else
        // Helpful for a preparation of the simulation or immediate execution (on cpu/gpu/locally, testing/running directly)
        custom_main();

#ifdef RUN_WITH_PYTHON_BACKEND
    mcmc::execution::finalize_python();
#endif
    return 0;
}

void custom_main()
{
    typedef lm_impl::lattice_system::ComplexXYModelParameters ModelParams;
    
    typedef std::complex<double> BasicType;
    typedef lm_impl::mcmc_update::ComplexLangevinUpdateParameters<ModelParams, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::ParallelUpdateWithAdpativeStepsizeParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "ComplexXYModelComplexLangevin";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {8, 8, 8};

    MCMCUpdateParams mcmc_update_parameters(0.001);

    double beta = 0.7;
    ModelParams model_parameters(beta, 0.0);

    // Thermalization time and measure interval in terms of Langevin time
    UpdateDynamicsParams update_dynamics_parameters(1.0, 0.01);

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Mean", "Config", "Energy"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    // Note that the measure interval is set to 1 since an update with adaptive step size is used
    // -> measure interval is determined in terms of Langein time in the UpdateDynamicsParams
    ExecutionParams execution_parameters(1, 20000, 0, {}, // optional additional measures
                                         {}, // Meausures which will be evaluated in terms of mean and error evaluation
                                         20); // Compute error based on Bootstrap method with 200 sampled sets of configurations 

    auto simulation_params = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, rel_data_path, "model_params", "mu", 0.0, std::sqrt(0.2), 2);
    simulation_params.write_to_file(rel_config_path);

    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);
}

// Rerun the simulation with ./ComplexXYModel expectation_value ComplexXYModelComplexLangevin
