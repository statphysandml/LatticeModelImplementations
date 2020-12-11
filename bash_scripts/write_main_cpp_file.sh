cat >"${src_path}/main.cpp" <<EOL
#ifndef MAIN_CPP
#define MAIN_CPP

EOL
if [ "$project_type" = "project" ]; then
cat >>"${src_path}/main.cpp" <<EOL
#include "../include/simulation_header.hpp"
#include "execution/executer.hpp"
EOL
else
cat >>"${src_path}/main.cpp" <<EOL
#include "simulation_header.hpp"
#include "execution/executer.hpp"
EOL
fi
cat >>"${src_path}/main.cpp" <<EOL

#endif

void custom_main();

int main(int argc, char **argv) {
    // initialize_executer_params(PROJECT_NAME, PYTHON_SCRIPTS_PATH, CLUSTER_MODE, CONDA_ACTIVATE_PATH, VIRTUAL_ENV); // optional

#ifdef PYTHON
    mcmc::execution::initialize_python();
#endif

    // A function of one of the first three if conditions is only called when an actual simulation takes place
    // (or for the generation of default parameters) based on a program that uses ./Main with arguments (from cpu/gpu/locally)
    if(argc > 1 and argc < 6)
    {
        run_from_file<typename from_file_simulation::MetropolisDynamics<XYModelParameters>::SystemBaseParams>(argc, argv);
    }
    else if(argc == 6)
    {
        from_file_simulation::run_based_on_algorithm<ComplexPolynomialModelParameters>(argc, argv);
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
}

void custom_main()
{
    typedef IsingModelParameters ModelParams;

    typedef double BasicType;
    typedef MetropolisUpdateParameters<ModelParams, IsingModelSampler> MCMCUpdateParams;
    typedef SequentialUpdateParameters UpdateDynamicsParams;
    typedef LatticeParameters<BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string target_name = "IsingModelMetropolis";
    std::string rel_config_path = "/configs/" + target_name + "/";
    std::string rel_data_path = "/data/" + target_name + "/";

    std::vector<double> dimensions {16, 16};

    MCMCUpdateParams mcmc_update_parameters(0.1);

    double beta = 0.4;
    ModelParams model_parameters(beta, 1.0, 0.0);

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Config", "AbsMean", "Mean", "Energy", "Std"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}},
            rel_config_path
    );

    typedef ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(1, 10000, 100, {}, {});

    auto simparams = SimulationParameters< SystemBaseParams , ExecutionParams >::generate_traceable_simulation(
            lattice_parameters, execution_parameters, rel_config_path, rel_data_path, "model_params", "beta", 0.1, 0.7, 10);

    // Execute the simulation
    execute< SystemBaseParams > (ExecutionParams::name(), target_name, "/./", true,
                                 Executer::local, true);

    /* typedef XYModelParameters ModelParams;

    typedef double BasicType;
    typedef MetropolisUpdateParameters<ModelParams, GaussianSampler> MCMCUpdateParams;
    typedef SequentialUpdateParameters UpdateDynamicsParams;
    typedef LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string target_name = "XYModelMetropolis";
    std::string rel_config_path = "/configs/" + target_name + "/";
    std::string rel_data_path = "/data/" + target_name + "/";

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
    execute< SystemBaseParams > (ExecutionParams::name(), target_name, "/./", true,
                                 Executer::local, false);
    */
}

// Rerun the simulation with ./$project_name expectation_value XYModelMetropolis

EOL
