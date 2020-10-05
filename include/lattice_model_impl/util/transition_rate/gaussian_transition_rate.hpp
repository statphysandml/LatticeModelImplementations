//
// Created by lukas on 28.09.20.
//

#ifndef COMPLEXMONTECARLO_GAUSSIAN_TRANSITION_RATE_HPP
#define COMPLEXMONTECARLO_GAUSSIAN_TRANSITION_RATE_HPP

#include "../transition_rate_base.hpp"

template<typename T, typename ModelCl, typename SamplerCl>
struct GaussianTransitionRate : TransitionRateBase<T, ModelCl, SamplerCl>
{
    GaussianTransitionRate<T, ModelCl, SamplerCl>(const ModelCl &model_, SamplerCl& sampler_, T state_= T(0.0, 0.0), double normalization_factor_ = 1.0, const std::string expansion_="full")
            : TransitionRateBase<T, ModelCl, SamplerCl>(model_, sampler_, state_, normalization_factor_, expansion_)
    {
        if(this->sampler.name() != required_sampler())
        {
            std::cout << "Sampler and transition rate do not coincide" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if(this->expansion == "first") {
            this->current_drift = this->model.get_drift_term(state_);
            function_ptr_on_expansion = &TransitionRateBase<T, ModelCl, SamplerCl>::get_first_order_action_diff;
        }
        else if(this->expansion == "second")
        {
            this->current_second_order_drift = this->model.get_second_order_drift_term(state_);
            function_ptr_on_expansion = &TransitionRateBase<T, ModelCl, SamplerCl>::get_second_order_action_diff;
        }
        else
            function_ptr_on_expansion = &TransitionRateBase<T, ModelCl, SamplerCl>::get_action_diff;
    }

    #ifdef THRUST
    __host__ __device__
    #endif
    double operator() (double x_real)
    {
        T x {x_real, this->state.imag()};

        auto action_diff = (this->*function_ptr_on_expansion)(x);

        double real_arg = -0.5 * std::pow(transformer(x.real()) - this->state.real(), 2.0) / (2.0 * this->eps) + action_diff.real();  // +
                          // 0.5 * std::pow(new_imag_state - state.imag(), 2.0) / (2.0 * eps);

        return jacobian(x_real) * exp(real_arg) / this->normalization_factor;
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_action_diff(T &proposed_site)
    {
        auto new_potential = this->model.get_potential({transformer(proposed_site.real()), proposed_site.imag()});
        return -0.5 * (new_potential - this->current_potential);
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_first_order_action_diff(T &proposed_site)
    {
        auto drift_term = this->model.get_drift_term({transformer(proposed_site.real()), proposed_site.imag()});
        return -0.5 * (proposed_site.real() - this->state.real()) * (drift_term - this->current_drift);
    }

#ifdef THRUST
    __host__ __device__
#endif
    T get_second_order_action_diff(T &proposed_site)
    {
        auto drift_term = this->model.get_drift_term({transformer(proposed_site.real()), proposed_site.imag()});
        auto second_order_drift_term = this->model.get_second_order_drift_term({transformer(proposed_site.real()), proposed_site.imag()});
        return -0.5 * (proposed_site.real() - this->state.real()) * (drift_term - this->current_drift) - 0.25 * std::pow(proposed_site.real() - this->state.real(), 2.0) * (second_order_drift_term - this->current_second_order_drift);
    }

    std::pair<double, double> get_imag_state_and_metropolis_acceptance_rate(T x)
    {
        x.imag(this->state.imag());
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());

        return {compute_new_imag_state(imag_arg, this->eps, x.real(), this->state), exp(-1.0 * (new_potential.real() - this->current_potential.real()))};
    }

    std::pair<double, double> get_imag_state_and_symmetric_acceptance_rate(T x)
    {
        x.imag(this->state.imag());
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());

        return {compute_new_imag_state(imag_arg, this->eps, x.real(), this->state),
                exp(-0.5 * std::pow(x.real() - this->state.real(), 2.0) / (2.0 * this->eps)  -
                    0.5 * (new_potential.real() - this->current_potential.real())) / this->normalization_factor};
        // return {compute_new_imag_state(imag_arg, eps, x.real(), state), exp(-1.0 * (new_potential.real() - current_potential.real()))};
    }

    #ifdef THRUST
    __host__ __device__
    #endif
    static double compute_new_imag_state(double imag_arg, const double eps, double x_real, const std::complex<double> &state)
    {
        if(x_real - state.real() == 0) // imag_arg is zero in this case
            return state.imag();
        else
            return state.imag() + 2.0 * eps * imag_arg / (x_real - state.real());
    }

    #ifdef THRUST
    __host__ __device__
    #endif
    double get_new_imag_state(T x)
    {
        x.imag(this->state.imag());
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());
        return compute_new_imag_state(imag_arg, this->eps, x.real(), this->state);
    }

    #ifdef THRUST
    __host__ __device__
    #endif
    double jacobian(const double x)
    {
        return -2.0 / (std::pow(x, 2.0) - 1.0);
    }

    struct transformer_func
    {
        #ifdef THRUST
        __host__ __device__
        #endif
        double operator() (const double val)
        {
            return std::log((1.0 + val) / (1.0 - val));
        }
    };

    static std::string required_sampler()
    {
        return "GaussianSampler";
    }

    transformer_func transformer;

    T (TransitionRateBase<T, ModelCl, SamplerCl>::*function_ptr_on_expansion)(T&);
};

#endif //COMPLEXMONTECARLO_GAUSSIAN_TRANSITION_RATE_HPP
