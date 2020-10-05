//
// Created by lukas on 30.09.20.
//

#ifndef COMPLEXMONTECARLO_COMPLEX_METROPOLIS_LANGEVIN_UPDATE_HPP
#define COMPLEXMONTECARLO_COMPLEX_METROPOLIS_LANGEVIN_UPDATE_HPP

#include <boost/math/tools/roots.hpp>

#include "../../mcmc_update_base.hpp"
#include "../../../util/integration/integration.h"

template<typename TransitionRate, typename ModelParameters, typename SamplerCl>
class ComplexMetropolisLangevinUpdate;

template<typename TransitionRate, typename ModelParameters, typename SamplerCl>
class ComplexMetropolisLangevinUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit ComplexMetropolisLangevinUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_)
    {}

    explicit ComplexMetropolisLangevinUpdateParameters(const double eps) : ComplexMetropolisLangevinUpdateParameters(json {{"eps", eps}})
    {}

    static std::string name() {
        return "ComplexMetropolisLangevinUpdate";
    }

    typedef ComplexMetropolisLangevinUpdate<TransitionRate, ModelParameters, SamplerCl> MCMCUpdate;

private:
    friend class ComplexMetropolisLangevinUpdate<TransitionRate, ModelParameters, SamplerCl>;
};


template<typename TransitionRate, typename ModelParameters, typename SamplerCl>
class ComplexMetropolisLangevinUpdate : public MCMCUpdateBase< ComplexMetropolisLangevinUpdate<TransitionRate, ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit ComplexMetropolisLangevinUpdate(
            const ComplexMetropolisLangevinUpdateParameters<TransitionRate, ModelParameters, SamplerCl> &up_,
            typename ModelParameters::Model & model_
    ) : MCMCUpdateBase< ComplexMetropolisLangevinUpdate<TransitionRate, ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_), model(model_), transition_rate(TransitionRate(model, this->sampler, std::complex<double>{0.0, 0.0}, 1.0, "full")), integrator(readdy::util::integration::integrator())
    {
        rand = std::uniform_real_distribution<double>(0.0, 1.0);
    }

    template<typename T>
    T operator() (const T site)
    {
        transition_rate.update_config(site);

        /* // Normalization
        auto normalization = compute_normalization_factor(transition_rate);
        transition_rate.set_normalization(normalization); */

        auto proposed_site = this->sampler.propose_state(site);

        auto result = transition_rate.get_imag_state_and_metropolis_acceptance_rate(proposed_site);

        auto new_imag_state = result.first;
        auto acceptance_factor = result.second;

        if(rand(gen) < std::min(1.0, acceptance_factor)) {
            return {proposed_site.real(), new_imag_state};
        }
        else {
            return {site.real(), new_imag_state};
        }
    }

private:
    const ComplexMetropolisLangevinUpdateParameters<TransitionRate, ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;
    TransitionRate transition_rate;

    readdy::util::integration::integrator integrator;
    std::uniform_real_distribution<double> rand;

    double compute_normalization_factor(TransitionRate &integrand)
    {
        auto normalization = integrand.compute_normalization_factor(integrator);
        if(normalization.second > 1e-10)
        {
            auto state = integrand.get_state();
            std::cout << "Large error estimation in normalization factor detected: " << normalization.second
                      << " for state " << state.real() << " + i" << state.imag() << std::endl;
        }
        return normalization.first;
    }
};

#endif //COMPLEXMONTECARLO_COMPLEX_METROPOLIS_LANGEVIN_UPDATE_HPP
