//
// Created by lukas on 28.02.20.
//

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

    typedef double BasicType;
    typedef lm_impl::lattice_system::IsingModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::lattice_system::IsingModelSampler> MCMCUpdateParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, MCMCUpdateParams, lm_impl::update_dynamics::SequentialUpdateParameters> SystemBaseParams;

    /* typedef double BasicType;
    typedef lm_impl::lattice_system::XYModelParameters ModelParams;
    typedef lm_impl::mcmc_update::HybridMonteCarloUpdateParameters<BasicType, ModelParams, mcmc::sampler::GaussianSampler> UpdateParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, UpdateParams, lm_impl::update_dynamics::GlobalLatticeUpdateParameters> SystemBaseParams; */

    /* typedef lm_impl::link::U1 BasicType;
    typedef lm_impl::link_lattice_system::UOneModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::link_lattice_system::UOneModelSampler> UpdateParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, UpdateParams, lm_impl::update_dynamics::SequentialUpdateParameters> SystemBaseParams; */

    /* typedef lm_impl::link::SU2<double> BasicType;
    typedef lm_impl::link_lattice_system::SUTwoModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::link_lattice_system::SUTwoModelSampler> UpdateParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, UpdateParams, lm_impl::update_dynamics::SequentialUpdateParameters> SystemBaseParams; */

    /* typedef lm_impl::link::ON<double, 4> BasicType;
    typedef lm_impl::lattice_system::ONModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::lattice_system::ONModelSampler> UpdateParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, UpdateParams, lm_impl::update_dynamics::SequentialUpdateParameters> SystemBaseParams; */
}

#endif //MAIN_SIMULATION_HEADER_HPP
