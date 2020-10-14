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
        auto integration_bounds = sampler.get_integration_bounds(state);
        lower_bound = integration_bounds.first;
        upper_bound = integration_bounds.second;

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

    virtual void update_config(T state_, double normalization_factor_ = 1.0)
    {
        state = state_;
        current_potential = model.get_potential(state_);
        normalization_factor = normalization_factor_;

        auto integration_bounds = sampler.get_integration_bounds(state);
        lower_bound = integration_bounds.first;
        upper_bound = integration_bounds.second;

        if(expansion == "first" or expansion == "second")
            current_drift = model.get_drift_term(state_);
        if(expansion == "second")
            current_second_order_drift = model.get_second_order_drift_term(state_);
    }

#ifdef THRUST
    __host__ __device__
#endif
    virtual double operator() (double x_real) = 0;

#ifdef THRUST
    __host__ __device__
#endif
    T get_action_diff(T &proposed_site)
    {
        auto new_potential = this->model.get_potential({transform(proposed_site.real()), proposed_site.imag()});
        return -0.5 * (new_potential - this->current_potential);
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_first_order_action_diff(T &proposed_site)
    {
        return -0.5 * (transform(proposed_site.real()) - state.real()) * current_drift;
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_second_order_action_diff(T &proposed_site)
    {
        return -0.5 * (transform(proposed_site.real()) - state.real()) * current_drift - 0.25 * std::pow(transform(proposed_site.real()) - state.real(), 2.0) * current_second_order_drift;
    }


    template<typename Integrator, typename Integrand>
    auto compute_normalization_factor(Integrator& integrator, Integrand& integrand)
    {
        return integrator.integrate(integrand, lower_bound, upper_bound, true);
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

    const double transform(const double val)
    {
        return sampler.transformer(val);
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