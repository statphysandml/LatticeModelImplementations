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
    typedef lm_impl::link::SU2<double> BasicType;
    typedef lm_impl::link_lattice_system::SUTwoModelParameters ModelParams;
    typedef lm_impl::mcmc_update::MetropolisUpdateParameters<ModelParams, lm_impl::link_lattice_system::SUTwoModelSampler> UpdateParams;
    typedef lm_impl::update_dynamics::SequentialUpdateParameters UpdateDynamicsParams;
    typedef lm_impl::lattice_system::LatticeParameters< BasicType, ModelParams, UpdateParams, UpdateDynamicsParams> SystemBaseParams;
}

#endif //MAIN_SIMULATION_HEADER_HPP

