//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_CUBIC_GAUSSIAN_MODEL_HPP
#define MAIN_CUBIC_GAUSSIAN_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


template<typename SamplerCl>
class CubicGaussianModel;

template<typename SamplerCl>
class CubicGaussianModelParameters : public SiteModelParameters {
public:
    explicit CubicGaussianModelParameters(const json params_) : SiteModelParameters(params_),
                                                     eps(get_value_by_key<double>("eps"))
    {}

    explicit CubicGaussianModelParameters(double eps_) : CubicGaussianModelParameters<SamplerCl>(json {
            {"eps", eps_}
    })
    {}

    static std::string name() {
        return "CubicGaussianModel";
    }

    typedef CubicGaussianModel<SamplerCl> Model;

private:
    friend class CubicGaussianModel<SamplerCl>;

    const double eps;
};

template<typename SamplerCl>
class CubicGaussianModel : public SiteModel<CubicGaussianModel<SamplerCl>>
{
public:
    explicit CubicGaussianModel(const CubicGaussianModelParameters<SamplerCl> &mp_) : mp(mp_), sampler(SamplerCl(mp.eps)) {}

    template<typename T>
    T random_state()
    {
        return sampler.template random_state<std::complex<double>>();
    }

    std::complex<double> propose_state(std::complex<double> site)
    {
        return sampler.template propose_state<std::complex<double>>(site);
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
    const CubicGaussianModelParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

// int CubicGaussianModel::init = 0;

#endif //MAIN_CUBIC_GAUSSIAN_MODEL_HPP
