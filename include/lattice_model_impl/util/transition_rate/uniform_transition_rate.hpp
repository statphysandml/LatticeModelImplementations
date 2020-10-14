//
// Created by lukas on 03.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_UNIFORM_TRANSITION_RATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UNIFORM_TRANSITION_RATE_HPP

#include "../transition_rate_base.hpp"

template<typename T, typename ModelCl, typename SamplerCl>
struct UniformTransitionRate : TransitionRateBase<T, ModelCl, SamplerCl>
{
    UniformTransitionRate<T, ModelCl, SamplerCl>(const ModelCl &model_, SamplerCl& sampler_, T state_= T(0.0, 0.0),
                                                     double normalization_factor_ = 1.0, const std::string expansion_="full")
            : TransitionRateBase<T, ModelCl, SamplerCl>(model_, sampler_, state_, normalization_factor_, expansion_)
    {}

    using TransitionRateBase<T, ModelCl, SamplerCl>::function_ptr_on_expansion;

#ifdef THRUST
    __host__ __device__
#endif
    double operator() (double x_real)
    {
        T x {x_real, this->state.imag()};

        auto action_diff = (this->*function_ptr_on_expansion)(x);

        double imag_arg = action_diff.imag();
        double real_arg = action_diff.real();

        double sampler_real_contribution = 1.0 / (2*this->eps) * (cos(imag_arg) + (pow(sin(imag_arg), 2.0) / cos(imag_arg)));

        if(sampler_real_contribution < 0.0) {
            std::cout << "Not real" << std::endl;
            // std::exit(EXIT_FAILURE);
        }
        return this->sampler.jacobian(x_real) * sampler_real_contribution * std::exp(real_arg) / this->normalization_factor;
    }

    std::pair<double, double> get_imag_state_and_metropolis_acceptance_rate(T x)
    {
        x.imag(this->state.imag());
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());

        T new_site_two = this->sampler.propose_state(this->state);
        double new_imag_state = compute_new_imag_state(imag_arg, this->eps, x.real(), new_site_two.real(), this->state);

        return {new_imag_state, exp(-1.0 * (new_potential.real() - this->current_potential.real()))};
    }

#ifdef THRUST
    __host__ __device__
#endif
    static double compute_new_imag_state(double imag_arg, const double eps, double x_real, double x_real_two, const std::complex<double> &state)
    {
        // state.imag() - (x_real - x_real_two) * tan(imag_arg_two);  // Both variants work - in this case imag_arg_two needs also to be calculated
        return state.imag() - (x_real_two - x_real) * tan(imag_arg);
    }

#ifdef THRUST
    __host__ __device__
#endif
    double get_new_imag_state(T x)
    {
        x.imag(this->state.imag());
        auto new_potential = this->model.get_potential(x);
        double imag_arg = -0.5 * (new_potential.imag() - this->current_potential.imag());

        T new_site_two = this->sampler.propose_state(this->state);
        return compute_new_imag_state(imag_arg, this->eps, x.real(), new_site_two.real(), this->state);
    }

    static std::string required_sampler()
    {
        return "UniformSampler";
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_UNIFORM_TRANSITION_RATE_HPP
