//
// Created by lukas on 05.10.20.
//

#ifndef COMPLEXMONTECARLO_ALGORITHMS_HPP
#define COMPLEXMONTECARLO_ALGORITHMS_HPP

#include "mcmc_simulation/header.hpp"

#include "../update_dynamics/update_dynamics_header.hpp"
#include "../mcmc_update/mcmc_update_header.hpp"

#include "../site/site_header.hpp"
#include "../lattice/lattice_header.hpp"
#include "../link_lattice/link_lattice_header.hpp"

#include "../util/integration_header.hpp"
#include "../util/implicit_integral_solver_header.hpp"
#include "../util/transition_rate_header.hpp"

namespace memory_site_algorithms
{
    struct MemorySiteAlgorithm
    {};

    template<typename ModelParameters>
    struct ComplexLangevinDynamics : MemorySiteAlgorithm
    {
        ComplexLangevinDynamics() = default;

        typedef std::complex<double> BasicType;
        typedef ComplexLangevinUpdateParameters<ModelParameters, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }

    };

    template<typename ModelParameters>
    struct SecondOrderComplexLangevinDynamics : MemorySiteAlgorithm
    {
        SecondOrderComplexLangevinDynamics() = default;

        typedef std::complex<double> BasicType;
        typedef ComplexLangevinSecondOrderUpdateParameters<ModelParameters, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }

    };

    template<typename ModelParameters>
    struct ComplexGaussianMetropolis
    {
        ComplexGaussianMetropolis() = default;

        typedef std::complex<double> BasicType;
        typedef GaussianTransitionRate<BasicType, typename ModelParameters::Model, mcmc::sampler::GaussianSampler> TransitionRate;
        typedef ComplexMetropolisLangevinUpdateParameters<TransitionRate, ModelParameters, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }
    };

    template<typename ModelParameters>
    struct ComplexHatFunctionMetropolis
    {
        ComplexHatFunctionMetropolis() = default;

        typedef std::complex<double> BasicType;
        typedef HatFunctionTransitionRate<BasicType, typename ModelParameters::Model, mcmc::sampler::HatFunctionSampler> TransitionRate;
        typedef ComplexMetropolisLangevinUpdateParameters<TransitionRate, ModelParameters, mcmc::sampler::HatFunctionSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }
    };

    template<typename ModelParameters, typename SamplerCl>
    struct ComplexUniformMetropolis
    {
        ComplexUniformMetropolis() = default;

        typedef std::complex<double> BasicType;
        typedef UniformTransitionRate<BasicType, typename ModelParameters::Model, SamplerCl> TransitionRate;
        typedef ComplexMetropolisLangevinUpdateParameters<TransitionRate, ModelParameters, SamplerCl> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }
    };

    template<typename ModelParameters>
    struct ComplexImlicitGaussianAlgorithm
    {
        ComplexImlicitGaussianAlgorithm() = default;

        typedef std::complex<double> BasicType;
        typedef GaussianTransitionRate<BasicType, typename ModelParameters::Model, mcmc::sampler::GaussianSampler> TransitionRate;
        typedef ComplexImplicitUpdateParameters<readdy::util::integration::integrator, iterative_solver, TransitionRate, ModelParameters, mcmc::sampler::GaussianSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }
    };

    template<typename ModelParameters>
    struct ComplexImlicitHatFunctionAlgorithm
    {
        ComplexImlicitHatFunctionAlgorithm() = default;

        typedef std::complex<double> BasicType;
        typedef HatFunctionTransitionRate<BasicType, typename ModelParameters::Model, mcmc::sampler::HatFunctionSampler> TransitionRate;
        typedef ComplexImplicitUpdateParameters<readdy::util::integration::integrator, iterative_solver, TransitionRate, ModelParameters, mcmc::sampler::HatFunctionSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }
    };

    template<typename ModelParameters, typename SamplerCl>
    struct ComplexImlicitUniformAlgorithm
    {
        ComplexImlicitUniformAlgorithm() = default;

        typedef std::complex<double> BasicType;
        typedef UniformTransitionRate<BasicType, typename ModelParameters::Model, SamplerCl> TransitionRate;
        typedef ComplexImplicitUpdateParameters<readdy::util::integration::integrator, iterative_solver, TransitionRate, ModelParameters, SamplerCl> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> SiteUpdateParameters;
        typedef SiteParameters<BasicType, ModelParameters, MCMCUpdateParams, SiteUpdateParameters> SystemBaseParams;

        static MCMCUpdateParams generate_update_parameters(const json params)
        {
            return MCMCUpdateParams(params);
        }
    };
}

#endif //COMPLEXMONTECARLO_ALGORITHMS_HPP
