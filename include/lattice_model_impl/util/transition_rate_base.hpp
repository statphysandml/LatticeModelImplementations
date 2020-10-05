//
// Created by lukas on 05.10.20.
//

#ifndef COMPLEXMONTECARLO_TRANSITION_RATE_BASE_HPP
#define COMPLEXMONTECARLO_TRANSITION_RATE_BASE_HPP

#include <complex>
#include <iostream>

template<typename T, typename ModelCl, typename SamplerCl>
struct TransitionRateBase
{
    TransitionRateBase<T, ModelCl, SamplerCl>(const ModelCl &model_, SamplerCl& sampler_, T state_= T(0.0, 0.0), double normalization_factor_ = 1.0, const std::string expansion_="full"):
            model(model_), sampler(sampler_), current_potential(model.get_potential(state_)), eps(sampler.get_eps()), state(state_), normalization_factor(normalization_factor_), expansion(expansion_)
    {
        if(expansion == "first") {
            current_drift = model.get_drift_term(state_);
            function_ptr_on_expansion = &TransitionRateBase<T, ModelCl, SamplerCl>::get_first_order_action_diff;
        }
        else if(expansion == "second")
        {
            current_second_order_drift = model.get_second_order_drift_term(state_);
            function_ptr_on_expansion = &TransitionRateBase<T, ModelCl, SamplerCl>::get_second_order_action_diff;
        }
        else
            function_ptr_on_expansion = &TransitionRateBase<T, ModelCl, SamplerCl>::get_action_diff;
    }

    void update_config(T state_, double normalization_factor_ = 1.0)
    {
        state = state_;
        current_potential = model.get_potential(state_);
        normalization_factor = normalization_factor_;

        if(expansion == "first" or expansion == "second")
            current_drift = model.get_drift_term(state_);
        if(expansion == "second")
            current_second_order_drift = model.get_second_order_drift_term(state_);
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_action_diff(T &proposed_site)
    {
        auto new_potential = model.get_potential({proposed_site.real(), proposed_site.imag()});
        return -0.5 * (new_potential - current_potential);
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_first_order_action_diff(T &proposed_site)
    {
        auto drift_term = model.get_drift_term({proposed_site.real(), proposed_site.imag()});
        return -0.5 * (proposed_site.real() - state.real()) * (drift_term - current_drift);
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_second_order_action_diff(T &proposed_site)
    {
        auto drift_term = model.get_drift_term({proposed_site.real(), proposed_site.imag()});
        auto second_order_drift_term = model.get_second_order_drift_term({proposed_site.real(), proposed_site.imag()});
        return -0.5 * (proposed_site.real() - state.real()) * (drift_term - current_drift) - 0.25 * std::pow(proposed_site.real() - state.real(), 2.0) * (second_order_drift_term - current_second_order_drift);
    }


    template<typename Integrator>
    auto compute_normalization_factor(Integrator& integrator)
    {
        return integrator.integrate(*this, lower_bound, upper_bound, true);
    }

    void set_normalization(const double normalization_factor_)
    {
        normalization_factor = normalization_factor_;
    }

    const T& get_state() const
    {
        return state;
    }

    const double get_lower_bound() const
    {
        return lower_bound;
    }

    const double get_upper_bound() const
    {
        return upper_bound;
    }

    const ModelCl &model;
    SamplerCl &sampler;

    T state;
    T current_potential;
    T current_drift;
    T current_second_order_drift;

    const double eps;
    double normalization_factor;

    const std::string expansion;
    double lower_bound;
    double upper_bound;

    T (TransitionRateBase<T, ModelCl, SamplerCl>::*function_ptr_on_expansion)(T&);
};

#endif //COMPLEXMONTECARLO_TRANSITION_RATE_BASE_HPP
