//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_SECOND_ORDER_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_SECOND_ORDER_UPDATE_HPP


#include "../../mcmc_update_base.hpp"


template<typename ModelParameters, typename SamplerCl>
class ComplexLangevinSecondOrderUpdate;


template<typename ModelParameters, typename SamplerCl>
class ComplexLangevinSecondOrderUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit ComplexLangevinSecondOrderUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                              epsilon(get_value_by_key<double>("epsilon")), sqrt2epsilon(sqrt(2 * get_value_by_key<double>("epsilon")))
    {}

    explicit ComplexLangevinSecondOrderUpdateParameters(
            const double epsilon_
    ) : ComplexLangevinSecondOrderUpdateParameters(json{
            {"epsilon", epsilon_},
            {"eps", epsilon_}
    })
    {}

    static std::string name() {
        return "ComplexLangevinSecondOrderUpdate";
    }

    typedef ComplexLangevinSecondOrderUpdate<ModelParameters, SamplerCl> MCMCUpdate;

private:
    friend class ComplexLangevinSecondOrderUpdate<ModelParameters, SamplerCl>;

    const double epsilon;
    const double sqrt2epsilon;
};


template<typename ModelParameters, typename SamplerCl>
class ComplexLangevinSecondOrderUpdate : public MCMCUpdateBase< ComplexLangevinSecondOrderUpdate<ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit ComplexLangevinSecondOrderUpdate(const ComplexLangevinSecondOrderUpdateParameters<ModelParameters, SamplerCl> &up_, typename ModelParameters::Model & model_)
        : MCMCUpdateBase< ComplexLangevinSecondOrderUpdate< ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_), model(model_)
    {
        normal = std::normal_distribution<double> (0,1);
    }

    template<typename T>
    T estimate_drift_term(const T site)
    {
        return model.get_drift_term(site);
    }

    template<typename T>
    T estimate_drift_term(const T site, const std::vector< T* > neighbours)
    {
        return model.get_drift_term(site, neighbours);
    }

    template<typename T>
    T operator() (const T site)
    {
        return update(site, model.get_drift_term(site), up.epsilon, up.sqrt2epsilon);
    }

    template<typename T>
    T operator() (const T site, const std::vector< T* > neighbours)
    {
        return update(site, model.get_drift_term(site, neighbours), up.epsilon, up.sqrt2epsilon);
    }

    template<typename T>
    T operator() (const T site, const double KMax, const double KExpectation)
    {
        T eps_drift_term = model.get_drift_term(site);
        double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
        return update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
    }

    template<typename T>
    T operator() (const T site, const std::vector< T* > neighbours, const double KMax, const double KExpectation)
    {
        T eps_drift_term = model.get_drift_term(site);
        double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
        return update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
    }

private:
    const ComplexLangevinSecondOrderUpdateParameters<ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;

    std::normal_distribution<double> normal;

    template<typename T>
    T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon)
    {
        T second_order_drift_term = model.get_second_order_drift_term(site);
        T new_site = {site.real(), site.imag()};
        new_site.real(site.real() - (epsilon * drift_term.real() + sqrt2epsilon * normal(gen)) / (1 + 0.5 * epsilon * second_order_drift_term.real()));
        new_site.imag(site.imag() - epsilon * drift_term.imag() - 0.5 * epsilon * (new_site.real() - site.real()) * second_order_drift_term.imag());

        // second_order_drift_term = model.get_second_order_drift_term(new_site); // Mixed
        // T drift_term_ = model.get_drift_term(new_site);
        // new_site.imag(site.imag() - epsilon * drift_term_.imag() - 0.5 * epsilon * (new_site.real() - site.real()) * second_order_drift_term.imag());

        // new_site.imag(site.imag() - epsilon * drift_term.imag() / (1 + 0.5 * epsilon * second_order_drift_term.imag())); // in z
        return new_site;
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_SECOND_ORDER_UPDATE_HPP
