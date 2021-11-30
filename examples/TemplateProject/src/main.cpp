#include "../include/TemplateProject/config.h"
#include "../include/TemplateProject/simulation_header.hpp"


void custom_main();

int main(int argc, char **argv) {
    param_helper::fs::prfs::set_relative_path_to_project_root_dir("../");

    // Initialization - Only needed for GPU and CPU runs
    mcmc::execution::initialize_executer_params(PROJECT_NAME, CLUSTER_MODE);

#ifdef RUN_WITH_PYTHON_BACKEND
    mcmc::execution::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    // A function of one of the first three if-conditions is only called when an actual simulation takes place (from cpu/gpu/locally)
    // (or for the generation of default parameters) based on a program that uses ./Main with arguments
    if(argc > 1 and argc < 6)
    {
        mcmc::execution::run_from_file< typename from_file_simulation::MetropolisDynamics<
            lm_impl::link::ON<double, 4>,
            lm_impl::lattice_system::ONModelSampler,
            lm_impl::lattice_system::ONModelParameters
        >::SystemBaseParams>(argc, argv);
    }
    else if(argc == 6)
    {
        from_file_simulation::run_based_on_algorithm< lm_impl::link::ON<double, 4>,
            lm_impl::lattice_system::ONModelSampler,
            lm_impl::lattice_system::ONModelParameters >(argc, argv);
    }
    else if(argc == 7)
    {
        from_file_simulation::run_based_on_model_and_algorithm(argc, argv);
    }
    else
        // Helpful for a preparation of the simulation or immediate execution (on cpu/gpu/locally, testing/running directly)
        custom_main();

#ifdef RUN_WITH_PYTHON_BACKEND
    mcmc::execution::finalize_python();
#endif
    return 0;
}

/* Example for running a simulation of the O(n) Model with a Metropolis algorithm -
 * On the example of the Execution mode: ExpectationValue
 * Further, the different ways to execute code are explained. */

void custom_main()
{
    typedef lm_impl::lattice_system::ONModelParameters ModelParams;

    typedef lm_impl::link::ON<double, 4> BasicType;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::lattice_system::ONModelSampler> MCMCUpdateParams;
    typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    std::string model_name = "ONModelMetropolis";
    std::string rel_config_path = "/configs/" + model_name + "/";
    std::string rel_data_path = "/data/" + model_name + "/";

    std::vector<int> dimensions {4, 4};

    MCMCUpdateParams mcmc_update_parameters(0.1);

    ModelParams model_parameters(json{
        {"kappa", 1.0},
        {"lambda", 1.0}
    });

    UpdateDynamicsParams update_dynamics_parameters;

    SystemBaseParams lattice_parameters(
            json {
                    {"dimensions", dimensions},
                    {"measures", {"Config", "Mean", "SecondMoment", "Energy"}},
                    {ModelParams::param_file_name(), model_parameters.get_json()},
                    {MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()},
                    {UpdateDynamicsParams::param_file_name(), update_dynamics_parameters.get_json()}}
    );

    typedef mcmc::execution::ExpectationValueParameters ExecutionParams;
    ExecutionParams execution_parameters(10, 10000, 10000, {}, // optional additional measures
                                         {"TwoPointCorrelation"}, // Meausures which will be evaluated in terms of mean and error evaluation
                                         200); // Compute error based on Bootstrap method with 200 sampled sets of configurations

    auto simulation_params = mcmc::simulation::SimulationParameters< SystemBaseParams , ExecutionParams >::generate_simulation(
            lattice_parameters, execution_parameters, rel_data_path, "model_params", "kappa", -1.0, 1.0, 21);

    // Store the simulation parameters - Only necessary if one wants to run the simulation again or on a CPU cluster, for example.
    simulation_params.write_to_file(rel_config_path);

    // Running the code with the executer function
    mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name);

    // Alternatively, you can set up and run the actual simulation. This also works without storing the simulation
    // parameters. Note, that in this case the expectation values are not computed in Python, as it is the case
    // running the simulation with the executer.
    /* mcmc::simulation::Simulation< SystemBaseParams, ExecutionParams > simulation(simulation_params);
     * simulation.run(); */

    // When parameters are stored, there are the following ways to run the simulation:

    // a)
    // Load the simulation parameters and run the simulation as before
    /* // Load simulation parameters from file - Note: It is important that the template parameters coincide for storing and loading
     * auto from_file_simulation_params = mcmc::simulation::SimulationParameters< SystemBaseParams, ExecutionParams >::generate_simulation_from_file(
            rel_config_path);
     * // Setting up and running the actual simulation
     * mcmc::simulation::Simulation< SystemBaseParams, ExecutionParams > from_file_simulation(from_file_simulation_params);
     * from_file_simulation.run(); */

    // b)
    // With the help of the executer. In this case, the expectation values are also computed in python if you have
    // python integrated - According to the ExpectationValue Mode - This is what is done in this program.
    /* mcmc::execution::execute< SystemBaseParams > (ExecutionParams::name(), model_name); */

    // c)
    // By running:
    //     ./TemplateProject expectation_value ONModelMetropolis
    // in the terminal. This only works if and only if you keep the respective if-function in the main function AND the
    // provided systembase template parameter coincides with the one in your config files. With the command exactly
    // the same as in b) will happen, with the exception that the SystemBaseParams parameter is passed in the main
    // function of the project/simulation.

    // d)
    // By changing further parameters of the execute command, you can submit a job to run the code on a cluster
    // (execute_code=true) or only prepare the respective bash script in the cpu_cluster_runs/ directory of the
    // project (execute_code=false). On the cluster, the simulation is executed based on the code in c). Accordingly, it is again
    // important that the restrictions to the correct system base parameters in c) hold. This behaviour can be changed by
    // changing the main function. To have higher flexibility, the execute function has a vector of strings as a last argument. This
    // allows passing additional parameters via the bash script to the main function of the project/simulation (of this file) when
    // the code is actually executed on the cluster.
    /* std::string sim_root_dir = "./";
     * bool rel_path = true;
     * bool execute_code = true; */  // Change to false to only generate the necessary bash_script for the execution on a cluster.
    // The bash_script is generated in the cpu_cluster_runs/ directory. In dependence on the cluster_mode variable,
    // the execute function either submits the job to the cluster or executes the code locally on your machine. The (terminal)
    // ouput of the simulation is in both cases copied to the cpu_cluster_runs/ directory.
    /* std::vector<std::string> additional_args = {};
     * mcmc::execution::execute< SystemBaseParams > (
     *        ExecutionParams::name(), model_name, sim_root_dir, rel_path,
     *        mcmc::execution::Executer::on_cpu_cluster, execute_code, additional_args); */

    // Additional remarks: By passing the cluster_mode with -DCLUSTER_MODE="local" or -DCLUSTER_MODE="on_cluster"
    // to cmake, you can switch between a "testing" mode, where the job is executed locally ("local") and an actual submission to
    // the cluster ("on_cluster"). For example, execute in the release directory:
    //     cmake ../cmake/ -DCMAKE_BUILD_TYPE=Release -DCLUSTER_MODE=on_cluster
    //     make -j4
}

// Rerun the simulation with ./TemplateProject expectation_value ONModelMetropolis
