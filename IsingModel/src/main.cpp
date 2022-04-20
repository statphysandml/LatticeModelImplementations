#include "../include/IsingModel/config.h"

#include <mcmc_simulation/header.hpp>
#include <mcmc_simulation/util/intervals.hpp>
#include <modes/mode_header.hpp>
#include <mcmc_simulation/util/random.hpp>

// #include <lattice_model_impl/representations/link_header.hpp>

#include <lattice_model_impl/update_dynamics/update_dynamics_header.hpp>
#include <lattice_model_impl/mcmc_method/mcmc_method_header.hpp>

// #include <lattice_model_impl/site/site_header.hpp>
#include <lattice_model_impl/lattice/lattice_header.hpp>
// #include <lattice_model_impl/link_lattice/link_lattice_header.hpp>


int main(int argc, char **argv) {
    // Initialize project dependent parameters
    param_helper::proj::set_relative_path_to_project_root_dir("../");

#ifdef PYTHON_BACKEND
    mcmc::util::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    //[ Defining the MCMC system and further important variables

    // Name of the simulation
    const std::string target_name = "IsingModelSimulation";

    // Directory for storing the results
    std::string rel_results_dir = "/results/" + target_name + "/";
    // Directory for storing the simulation data
    std::string rel_data_dir = "/data/" + target_name + "/";

    typedef lm_impl::lattice_system::IsingModel MCMCModel;
    
    typedef double BasicType;
    typedef lm_impl::lattice_system::IsingModelSampler Sampler;
    typedef lm_impl::mcmc_method::MetropolisUpdate<MCMCModel, Sampler> MCMCMethod;
    typedef lm_impl::update_dynamics::SequentialUpdate UpdateDynamics;
    typedef lm_impl::lattice_system::LatticeSystem<BasicType, MCMCModel, MCMCMethod, UpdateDynamics, Sampler> Lattice;

    std::vector<int> dimensions{4, 4};

    Sampler sampler;

    MCMCModel model(0.4, 1.0, 0.0);

    MCMCMethod mcmc_method;

    UpdateDynamics update_dynamics;

    // Setting up the system
    Lattice system(
        sampler,
        model,
        mcmc_method,
        update_dynamics,
        dimensions,
        "nearest_neighbour"
    );

    // Setting up measurement processor
    typedef mcmc::measures::ReadableMeasure ReadableMeasureProcessor;
    ReadableMeasureProcessor readable_measures(rel_data_dir);

    //]


    //[ Equilibrium time simulation

    // Setting up equilibrium time parameters
    typedef mcmc::mode::EquilibriumTime EquilibriumTimeParams;
    EquilibriumTimeParams equilibrium_time_parameters(
        10, // sample_size
        100, // number_of_steps
        0.005, // confidence_range
        10, // confidence_window
        "Mean" // measure
    );

    // Prepare the simulation
    auto beta_intervals = mcmc::util::linspace(0.1, 0.7, 5);
    auto equilibrium_time_simulation = mcmc::simulation::Simulation<
        Lattice, EquilibriumTimeParams, ReadableMeasureProcessor>::generate_simulation(
            system,
            equilibrium_time_parameters,
            readable_measures,
            "mcmc_model", // running_parameter_kind
            "beta", // running parameter (rp)
            beta_intervals // rp_intervals
    );
    
    // Run and evaluate the simulation
    equilibrium_time_simulation.run();
    equilibrium_time_simulation.eval(rel_results_dir);

    //]

    
    //[ Correlation time simulation

    // Relative path to the previously computed equilibrium time results
    std::string rel_equilibrium_time_results_path = "/results/" + target_name + "/";

    // Setting up correlation time parameters
    typedef mcmc::mode::CorrelationTime CorrelationTimeParams;
    CorrelationTimeParams correlation_time_parameters(
        1000, // minimum_sample_size
        400, // maximum_correlation_time
        rel_equilibrium_time_results_path, // equilibrium_time_rel_results_dir
        "Mean", // measure
        "hot" // starting_mode
    );
    
    // Prepare the simulation
    auto correlation_time_simulation = mcmc::simulation::Simulation<
        Lattice, CorrelationTimeParams, ReadableMeasureProcessor>::generate_simulation(
            system,
            correlation_time_parameters,
            readable_measures,
            "mcmc_model", // running_parameter_kind
            "beta", // running parameter (rp)
            beta_intervals // rp_intervals
    );

    // Run and evaluate the simulation
    correlation_time_simulation.run();
    correlation_time_simulation.eval(rel_results_dir);

    //]


    //[ Expectation Value

    // Relative path to the previously computed correlation time results
    std::string rel_correlation_time_results_path = "/results/" + target_name + "/";

    // Setting up expectation value parameters
    typedef mcmc::mode::ExpectationValue ExpectationValueParams;
    ExpectationValueParams expectation_value_parameters(
        rel_correlation_time_results_path, // correlation_time_rel_results_dir
        1000, //  number_of_measurements
        rel_equilibrium_time_results_path, // equilibrium_time_rel_results_dir
        {"Config", "Mean", "AbsMean", "SecondMoment", "Action", "AcceptanceRate", "EnergyViolation",
         "ExponentialEnergyViolation"}, // measures
         {}, // post_measures
         "hot", // starting_mode
         "statistical" // error_type
    );

    // Prepare the simulation
    auto expectation_value_simulation = mcmc::simulation::Simulation<
            Lattice, ExpectationValueParams, ReadableMeasureProcessor>::generate_simulation(
            system,
            expectation_value_parameters,
            readable_measures,
            "mcmc_model", // running_parameter_kind
            "beta", // running parameter (rp)
            beta_intervals // rp_intervals
    );
    
    // Run and evaluate the simulation
    expectation_value_simulation.run();
    expectation_value_simulation.eval(rel_results_dir);

    //]

    // Finalization
#ifdef PYTHON_BACKEND
    mcmc::util::finalize_python();
#endif
    return 0;
}