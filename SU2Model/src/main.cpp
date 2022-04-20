#include "../include/SU2Model/config.h"

#include <lattice_model_impl/representations/link_header.hpp>

#include <mcmc_simulation/header.hpp>
#include <mcmc_simulation/util/intervals.hpp>
#include <modes/mode_header.hpp>
#include <mcmc_simulation/util/random.hpp>

#include <lattice_model_impl/update_dynamics/update_dynamics_header.hpp>
#include <lattice_model_impl/mcmc_method/mcmc_method_header.hpp>

// #include <lattice_model_impl/site/site_header.hpp>
#include <lattice_model_impl/lattice/lattice_header.hpp>
#include <lattice_model_impl/link_lattice/link_lattice_header.hpp>

#include <lattice_model_impl/sampler/link_sampler.hpp>


int main(int argc, char **argv) {
    // Initialize project dependent parameters
    param_helper::proj::set_relative_path_to_project_root_dir("../");

#ifdef PYTHON_BACKEND
    mcmc::util::initialize_python(PYTHON_SCRIPTS_PATH);
#endif

    //[ Defining the MCMC system and further important variables

    // Name of the simulation
    const std::string target_name = "SU2ModelMetropolis";

    // Directory for storing the results
    std::string rel_results_dir = "/results/" + target_name + "/";
    // Directory for storing the simulation data
    std::string rel_data_dir = "/data/" + target_name + "/";

    typedef lm_impl::link_lattice_system::SU2Model MCMCModel;

    typedef lm_impl::link::SU2<double> BasicType;
    typedef lm_impl::link_lattice_system::LinkSampler Sampler;
    typedef lm_impl::mcmc_method::MetropolisUpdate<MCMCModel, Sampler> MCMCMethod;
    typedef lm_impl::update_dynamics::SequentialUpdate UpdateDynamics;
    typedef lm_impl::lattice_system::LatticeSystem<BasicType, MCMCModel, MCMCMethod, UpdateDynamics, Sampler> Lattice;

    std::vector<int> dimensions {4, 4, 4, 4};
    
    Sampler sampler(0.45);

    MCMCModel model(2.3);

    MCMCMethod mcmc_method;

    UpdateDynamics update_dynamics;

    // Setting up the system
    Lattice system(
        sampler,
        model,
        mcmc_method,
        update_dynamics,
        dimensions,
        "plaquette"
    );

    // Setting up measurement processor
    typedef mcmc::measures::ReadableMeasure ReadableMeasureProcessor;
    ReadableMeasureProcessor readable_measures(rel_data_dir);

    //]

    //[ Expectation Value

    typedef mcmc::mode::ExpectationValue ExpectationValueParams;
    ExpectationValueParams expectation_value_parameters(
        20, // correlation_time_rel_results_dir
        500, //  number_of_measurements
        100, // equilibrium_time_rel_results_dir
        {"Config", "AveragePlaquetteAction"}, // measures
        {}, // post_measures
        "hot", // starting_mode
        "statistical" // error_type
    );

    // Prepare the simulation
    auto beta_intervals = mcmc::util::linspace(1.7, 2.9, 5);
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

    //] */

    // Finalization
#ifdef PYTHON_BACKEND
    mcmc::util::finalize_python();
#endif
    return 0;
}
