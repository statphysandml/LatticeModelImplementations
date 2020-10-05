//
// Created by lukas on 03.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_UNIFORM_TRANSITION_RATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UNIFORM_TRANSITION_RATE_HPP

#include "../transition_rate_base.hpp"

template<typename T, typename ModelCl, typename SamplerCl>
struct UniformTransitionRate : TransitionRateBase<T, ModelCl, SamplerCl>
{
    using TransitionRateBase<T, ModelCl, SamplerCl>::TransitionRateBase;
    using TransitionRateBase<T, ModelCl, SamplerCl>::function_ptr_on_expansion;

#ifdef THRUST
    __host__ __device__
#endif
    double operator() (double x_real)
    {
        T x {x_real, this->state.imag()};

        double sign = get_sign(x);

        auto action_diff = (this->*function_ptr_on_expansion)(x);

        double imag_arg = action_diff.imag();
        double real_arg = action_diff.real();

        double new_imag_state = compute_new_imag_state(imag_arg, this->eps, x.real(), this->state, sign);

        double sampler_real_contribution = 1.0/this->eps * (1.0 - sign * (x.real() - this->state.real()) / this->eps) * cos(imag_arg) +
                                           1.0/this->eps * sign * (new_imag_state - this->state.imag()) / this->eps * sin(imag_arg);

        /* double sampler_imag_contribution = 1.0/eps * (1.0 - sign * (x.real() - state.real()) / eps) * sin(imag_arg) -
                                           1.0/eps * sign * (new_imag_state - state.imag()) / eps * cos(imag_arg); */

        if(sampler_real_contribution < 0.0) {
            std::cout << "Not real" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return sampler_real_contribution * std::exp(real_arg) / this->normalization_factor;
    }

    std::pair<double, double> get_imag_state_and_metropolis_acceptance_rate(T x)
    {
        x.imag(this->state.imag());
        double sign = get_sign(x);
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());

        double new_imag_state = compute_new_imag_state(imag_arg, this->eps, x.real(), x.real(), this->state);

        return {new_imag_state, exp(-1.0 * (new_potential.real() - this->current_potential.real()))};
    }

#ifdef THRUST
    __host__ __device__
#endif
    static double compute_new_imag_state(double imag_arg, const double eps, double x_real, double x_real_two, const std::complex<double> &state)
    {
        return state.imag() - (x_real - x_real_two) * tan(imag_arg); // -0.5 * (new_potential_two.imag() - current_potential.imag())
    }

#ifdef THRUST
    __host__ __device__
#endif
    double get_new_imag_state(T x)
    {
        x.imag(this->state.imag());
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());

        return compute_new_imag_state(imag_arg, this->eps, x.real(), x.real(), this->state);
    }

#ifdef THRUST
    __host__ __device__
#endif
    double get_sign(T x)
    {
        if(x.real() - this->state.real() == 0)
            return 1.0;
        else
            return (x.real() - this->state.real()) / std::abs(x.real() - this->state.real());
    }

    struct transformer_func
    {
#ifdef THRUST
        __host__ __device__
#endif
        double operator() (const double val)
        {
            return val;
        }
    };

    static std::string required_sampler()
    {
        return "UniformSampler";
    }

    transformer_func transformer;
};

#endif //LATTICEMODELIMPLEMENTATIONS_UNIFORM_TRANSITION_RATE_HPP
