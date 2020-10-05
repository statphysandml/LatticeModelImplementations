//
// Created by lukas on 01.10.20.
//

#ifndef COMPLEXMONTECARLO_COMPLEX_UNIFORM_UPDATE_HPP
#define COMPLEXMONTECARLO_COMPLEX_UNIFORM_UPDATE_HPP


#include "../../mcmc_update_base.hpp"


template<typename ModelParameters, typename SamplerCl>
class ComplexUniformUpdate;


template<typename ModelParameters, typename SamplerCl>
class ComplexUniformUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit ComplexUniformUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                   epsilon(get_value_by_key<double>("epsilon"))
    {}

    explicit ComplexUniformUpdateParameters(
            const double epsilon_
    ) : ComplexUniformUpdateParameters(json {
            {"epsilon", epsilon_},
            {"eps", epsilon_}
    })
    {}

    static std::string name() {
        return "ComplexUniformUpdate";
    }

    typedef ComplexUniformUpdate<ModelParameters, SamplerCl> MCMCUpdate;

private:
    friend class ComplexUniformUpdate<ModelParameters, SamplerCl>;

    const double epsilon;
};


template<typename ModelParameters, typename SamplerCl>
class ComplexUniformUpdate : public MCMCUpdateBase< ComplexUniformUpdate<ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit ComplexUniformUpdate(const ComplexUniformUpdateParameters<ModelParameters, SamplerCl> &up_, typename ModelParameters::Model &model_)
        : MCMCUpdateBase<ComplexUniformUpdate<ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_), model(model_)
    {
        uniform = std::uniform_real_distribution<double> (-up.epsilon, up.epsilon);
        rand = std::uniform_real_distribution<double> (0.0, 1.0);
    }

    template<typename T>
    T operator() (const T site)
    {
        T new_site = this->sampler.propose_state(site);

        auto current_potential = model.get_potential(site);
        auto new_potential = model.get_potential(new_site);

        T new_site_two = this->sampler.propose_state(site);

        auto new_imag_state = site.imag() - (new_site_two.real() - new_site.real()) * tan(-0.5 * (new_potential.imag() - current_potential.imag()));

        /* auto new_potential_two = model.get_potential(new_site_two); // Both variants work
        auto new_imag_state = site.imag() - (new_site.real() - new_site_two.real()) * tan(-0.5 * (new_potential_two.imag() - current_potential.imag())); */

        if(rand(gen) < std::min(1.0, exp(-1.0 * (new_potential.real() - current_potential.real()))))
        {
            return {new_site.real(), new_imag_state};
        }
        else
        {
            return {site.real(), new_imag_state};
        }
    }

private:
    const ComplexUniformUpdateParameters<ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;

    std::uniform_real_distribution<double> uniform;
    std::uniform_real_distribution<double> rand;
};

#endif //COMPLEXMONTECARLO_COMPLEX_UNIFORM_UPDATE_HPP
