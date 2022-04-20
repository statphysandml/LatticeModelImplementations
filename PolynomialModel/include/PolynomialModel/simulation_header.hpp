#ifndef MAIN_SIMULATION_HEADER_HPP
#define MAIN_SIMULATION_HEADER_HPP

#ifdef RUN_WITH_PYTHON_BACKEND
#include <Python.h>
#endif

#include <lattice_model_impl/representations/link_header.hpp>

#include <mcmc_simulation/header.hpp>
#include <execution/executer.hpp>

#include <lattice_model_impl/update_dynamics/update_dynamics_header.hpp>
#include <lattice_model_impl/mcmc_method/mcmc_method_header.hpp>

#include <lattice_model_impl/site/site_header.hpp>
#include <lattice_model_impl/lattice/lattice_header.hpp>
#include <lattice_model_impl/link_lattice/link_lattice_header.hpp>


namespace from_file_simulation {
    typedef double BasicType;
    typedef lm_impl::site_system::PolynomialModelParameters ModelParams;
    typedef lm_impl::mcmc_method::LangevinUpdateParameters<ModelParams, lm_impl::sampler::GaussianSampler> MCMCMethodParams;
    typedef lm_impl::site_system::SiteParameters< BasicType, ModelParams, MCMCMethodParams> SystemBaseParams;
}

#endif //MAIN_SIMULATION_HEADER_HPP

