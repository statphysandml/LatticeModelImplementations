cat >"${src_path}/main.cpp" <<EOL
#ifndef MAIN_CPP
#define MAIN_CPP

EOL
if [ "$project_type" = "project" ]; then
cat >>"${src_path}/main.cpp" <<EOL
#include "../include/config.h"
#include "../include/simulation_header.hpp"
EOL
else
cat >>"${src_path}/main.cpp" <<EOL
#include "config.h"
#include "simulation_header.hpp"
EOL
fi
cat >>"${src_path}/main.cpp" <<EOL

#endif

void custom_main();

int main(int argc, char **argv) {
    param_helper::fs::prfs::set_relative_path_to_project_root_dir("../");

    // Initialization - Only needed for GPU and CPU runs
    mcmc::execution::initialize_executer_params(PROJECT_NAME, CLUSTER_MODE);

#ifdef PYTHON
    mcmc::execution::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    // A function of one of the first three if-conditions is only called when an actual simulation takes place
    // (or for the generation of default parameters) based on a program that uses ./Main with arguments (from cpu/gpu/locally)
    if(argc > 1 and argc < 6)
    {
        mcmc::execution::run_from_file<typename from_file_simulation::MetropolisDynamics<lm_impl::lattice_system::XYModelParameters>::SystemBaseParams>(argc, argv);
    }
    else if(argc == 6)
    {
        from_file_simulation::run_based_on_algorithm<lm_impl::site_system::ComplexPolynomialModelParameters>(argc, argv);
    }
    else if(argc == 7)
    {
        from_file_simulation::run_based_on_model_and_algorithm(argc, argv);
    }
    else
        // Helpful for a preparation of the simulation or immediate execution (on cpu/gpu/locally, testing/running directly)
        custom_main();

#ifdef PYTHON
    mcmc::execution::finalize_python();
#endif
    return 0;
}

void custom_main()
{
    typedef lm_impl::lattice_system::IsingModelParameters ModelParams;

    typedef double BasicType;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::lattice_system::IsingModelSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "IsingModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.1);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Config", "Mean", "Std"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 10000, 10000, {}, // optional additional measures
                                         {"AbsMean", "Energy"}); // Meausures which will be evaluated in terms of mean and error evaluation

    auto simparams = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, rel_data_path, "model_params", "beta", 0.1, 0.625, 21);
    simparams.write_to_file(rel_config_path);

    // Execute the simulation
    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name, "./", true, mcmc::execution::Executer::local, true);

    /* typedef XYModelParameters ModelParams;

    typedef double BasicType;
    typedef MetropolisUpdateParameters<ModelParams, GaussianSampler> MCMCUpdateParams;
    typedef SequentialUpdateParameters UpdateDynamicsParams;
    typedef LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "XYModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<double> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.1);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0, 0.1);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"XYMagnetization"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(100, 10000, 100, {}, {});

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            lattice_parameters, execution_parameters, rel_config_path, rel_data_path, "model_params", "beta", 0.05, 1.55, 10);

    // Execute the simulation
    execute< SystemBaseParams > (ExecutionParams::name(), model_name, "/./", true,
                                 Executer::local, false);
    */
}

// Rerun the simulation with ./$project_name expectation_value XYModelMetropolis

EOL
