//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_COMPLEX_CUBIC_MODEL_HPP
#define MAIN_COMPLEX_CUBIC_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


template<typename SamplerCl>
class ComplexCubicModel;

template<typename SamplerCl>
class ComplexCubicModelParameters : public SiteModelParameters {
public:
    explicit ComplexCubicModelParameters(const json params_) : SiteModelParameters(params_),
                                                                eps(get_value_by_key<double>("eps", 0.1))
    {}

    explicit ComplexCubicModelParameters(double eps_) : ComplexCubicModelParameters<SamplerCl>(json {
            {"eps", eps_}
    })
    {}

    static std::string name() {
        return "ComplexCubicModel";
    }

    typedef ComplexCubicModel<SamplerCl> Model;

private:
    friend class ComplexCubicModel<SamplerCl>;

    const double eps;
};

template<typename SamplerCl>
class ComplexCubicModel : public SiteModel<ComplexCubicModel<SamplerCl>>
{
public:
    explicit ComplexCubicModel(const ComplexCubicModelParameters<SamplerCl> &mp_) : mp(mp_), sampler(SamplerCl(mp.eps)) {}

    template<typename T>
    T random_state()
    {
        return sampler.template random_state<T>();
    }

    template<typename T>
    T propose_state(T site)
    {
        return sampler.template propose_state<T>(site);
    }

    /* std::complex<double> propose_state(const std::complex<double> site, const double KMax, const double KExpectation)
    {
        double eps = std::min(mp.eps, mp.eps * KExpectation / KMax);
        std::complex<double> state = site + eps * std::complex<double>(propose_normal(gen), propose_normal(gen));
        return state;
    } */

    static std::complex<double> get_drift_term(const std::complex<double> site)
    {
        return  {-2.0 * site.real() * site.imag(), -1.0 * (std::pow(site.imag(), 2) - std::pow(site.real(), 2))};
    }

    static  std::complex<double> get_second_order_drift_term(const std::complex<double> site)
    {
        return {-2.0 * site.imag(), 2.0 * site.real()};
    }

    static std::complex<double> get_potential(const std::complex<double> site)
    {
        // return std::complex<double>{0, 1} * std::pow(site, 3) / 3.0;
        return {-1.0 * std::pow(site.real(), 2) * site.imag() + std::pow(site.imag(), 3) / 3.0, std::pow(site.real(), 3) / 3.0 - std::pow(site.imag(), 2) * site.real() };
    }

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

private:
    const ComplexCubicModelParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

// int ComplexCubicModel::init = 0;

#endif //MAIN_COMPLEX_CUBIC_MODEL_HPP