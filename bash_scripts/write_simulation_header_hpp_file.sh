cat >"${include_path}/simulation_header.hpp" <<EOL
#ifndef MAIN_SIMULATION_HEADER_HPP
#define MAIN_SIMULATION_HEADER_HPP

#ifdef PYTHON
#include <Python.h>
#endif

#include "lattice_model_impl/representations/link_header.hpp"

#include "mcmc_simulation/header.hpp"
#include "execution/executer.hpp"

#include "lattice_model_impl/update_dynamics/update_dynamics_header.hpp"
#include "lattice_model_impl/mcmc_update/mcmc_update_header.hpp"

#include "lattice_model_impl/site/site_header.hpp"
#include "lattice_model_impl/lattice/lattice_header.hpp"
#include "lattice_model_impl/link_lattice/link_lattice_header.hpp"


namespace from_file_simulation {
    /* typedef double BasicType;
    typedef lm_impl::lattice_system::XYModelParameters ModelParams;
    typedef lm_impl::mcmc_update::HybridMonteCarloUpdateParameters<BasicType, ModelParams, mcmc::sampler::GaussianSampler> UpdateParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, UpdateParams, lm_impl::update_dynamics::GlobalLatticeUpdateParameters> SystemBaseParams; */

    template<typename ModelParameters>
    struct ComplexLangevinDynamics
    {
        typedef std::complex<double> BasicType;
        typedef lm_impl::mcmc_update::ComplexLangevinUpdateParameters<ModelParameters, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
        typedef lm_impl::update_dynamics::MemorySiteSimpleUpdateParameters<BasicType> UpdateDynamicsParams;
        typedef lm_impl::site_system::SiteParameters< BasicType, ModelParameters, MCMCUpdateParams, UpdateDynamicsParams > SystemBaseParams;
    };

    template<typename ModelParameters>
    struct MetropolisDynamics
    {
        typedef double BasicType;
        typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParameters, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
        typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
        typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParameters, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    };

    template<typename ModelParameters>
    void run_based_on_algorithm(int argc, char **argv)
    {
        std::string algorithm = std::string(argv[5]);

        if(algorithm == "ComplexLangevinDynamics")
            mcmc::execution::run_from_file<typename ComplexLangevinDynamics<ModelParameters>::SystemBaseParams>(argc, argv);
    }

    void run_based_on_model_and_algorithm(int argc, char **argv)
    {
        std::string model = std::string(argv[6]);

        if(model == "ComplexPolynomialModel")
            run_based_on_algorithm<lm_impl::site_system::ComplexPolynomialModelParameters>(argc, argv);
    }
}

#endif //MAIN_SIMULATION_HEADER_HPP

EOL