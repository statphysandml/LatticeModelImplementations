#ifndef MAIN_SIMULATION_HEADER_HPP
#define MAIN_SIMULATION_HEADER_HPP

#ifdef RUN_WITH_PYTHON_BACKEND
#include <Python.h>
#endif

#include <lattice_model_impl/representations/link_header.hpp>

#include <mcmc_simulation/header.hpp>
#include <execution/executer.hpp>

#include <lattice_model_impl/update_dynamics/update_dynamics_header.hpp>
#include <lattice_model_impl/mcmc_update/mcmc_update_header.hpp>

#include <lattice_model_impl/site/site_header.hpp>
#include <lattice_model_impl/lattice/lattice_header.hpp>
#include <lattice_model_impl/link_lattice/link_lattice_header.hpp>



namespace from_file_simulation {
    template<typename BasicType, typename Sampler, typename ModelParameters>
    struct MetropolisDynamics
    {
        typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParameters, Sampler> MCMCUpdateParams;
        typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
        typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParameters, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    };

    template<typename BasicType, typename Sampler, typename ModelParameters>
    void run_based_on_algorithm(int argc, char **argv)
    {
        std::string algorithm = std::string(argv[5]);

        if(algorithm == "MetropolisDynamics")
            mcmc::execution::run_from_file<typename MetropolisDynamics<BasicType, Sampler, ModelParameters>::SystemBaseParams>(argc, argv);
    }

    void run_based_on_model_and_algorithm(int argc, char **argv)
    {
        std::string model = std::string(argv[6]);

        if(model == "ONModel")
            run_based_on_algorithm< lm_impl::link::ON<double, 4>,
                lm_impl::lattice_system::ONModelSampler,
                lm_impl::lattice_system::ONModelParameters >(argc, argv);
    }
}

#endif //MAIN_SIMULATION_HEADER_HPP

