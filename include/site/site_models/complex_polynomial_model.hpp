//
// Created by lukas on 02.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


template<typename SamplerCl>
class ComplexPolynomialModel;

template<typename SamplerCl>
class ComplexPolynomialModelParameters : public SiteModelParameters {
public:
    explicit ComplexPolynomialModelParameters<SamplerCl>(const json params_) : SiteModelParameters(params_),
                                                                               eps(get_value_by_key<double>("eps", 0.1)),
                                                                               sigma(get_value_by_key< std::complex<double> >("sigma")),
                                                                               lambda(get_value_by_key< std::complex<double> >("lambda")),
                                                                               h(get_value_by_key< std::complex<double> >("h", 0.0))
    {}

    explicit ComplexPolynomialModelParameters<SamplerCl>(std::complex<double> lambda_, std::complex<double> sigma_, std::complex<double> h_, double eps_) : ComplexPolynomialModelParameters(json {
            {"lambda", lambda_},
            {"sigma", sigma_},
            {"h", h_},
            {"eps", eps_}
    })
    {}

    static std::string name() {
        return "ComplexPolynomialModel";
    }

    typedef ComplexPolynomialModel<SamplerCl> Model;

private:
    friend class ComplexPolynomialModel<SamplerCl>;

    const double eps;
    const std::complex<double> lambda;
    const std::complex<double> sigma;
    const std::complex<double> h;
};

template<typename SamplerCl>
class ComplexPolynomialModel : public SiteModel< ComplexPolynomialModel<SamplerCl> >
{
public:
    explicit ComplexPolynomialModel(const ComplexPolynomialModelParameters<SamplerCl> &mp_) : mp(mp_), sampler(SamplerCl(mp.eps)) {}

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

    std::complex<double> get_potential(const std::complex<double> site) const
    {
        return 0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site;
    }

    double get_potential(double site) const
    {
        return (0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site).real();
    }

    double get_imag_potential(double site) const
    {
        return (0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site).imag();
    }

    /* template<typename func>
    std::complex<double> get_potential(const std::complex<double> site, func& transformer) const
    {
        std::complex<double> site_{transformer(site.real()), site.imag()};
        return 0.5 * mp.sigma * std::pow(site_, 2) + 0.25 * mp.lambda * std::pow(site_, 4) + mp.h * site_;
    } */

    std::complex<double> get_drift_term(const std::complex<double> site) const
    {
        return mp.sigma * site + mp.lambda * std::pow(site, 3) + mp.h;
        /* return  {mp.sigma.real() * site.real() - mp.sigma.imag() * site.imag() +
                    mp.lambda.real() * std::pow(site.real(), 3) - 3.0 * mp.lambda.real() * std::pow(site.imag(), 2) * site.real() +
                    mp.lambda.imag() * std::pow(site.imag(), 3) - 3.0 * mp.lambda.imag() * std::pow(site.real(), 2) * site.imag(),
                 mp.sigma.imag() * site.real() + mp.sigma.real() * site.imag() -
                    mp.lambda.real() * std::pow(site.imag(), 3) + 3.0 * mp.lambda.real() * std::pow(site.real(), 2) * site.imag() +
                    mp.lambda.imag() * std::pow(site.real(), 3) - 3.0 * mp.lambda.imag() * std::pow(site.imag(), 2) * site.real()}; */
    }

    double get_drift_term(const double site) const
    {
        return (mp.sigma * site + mp.lambda * std::pow(site, 3) + mp.h).imag();
        /* return  {mp.sigma.real() * site.real() - mp.sigma.imag() * site.imag() +
                    mp.lambda.real() * std::pow(site.real(), 3) - 3.0 * mp.lambda.real() * std::pow(site.imag(), 2) * site.real() +
                    mp.lambda.imag() * std::pow(site.imag(), 3) - 3.0 * mp.lambda.imag() * std::pow(site.real(), 2) * site.imag(),
                 mp.sigma.imag() * site.real() + mp.sigma.real() * site.imag() -
                    mp.lambda.real() * std::pow(site.imag(), 3) + 3.0 * mp.lambda.real() * std::pow(site.real(), 2) * site.imag() +
                    mp.lambda.imag() * std::pow(site.real(), 3) - 3.0 * mp.lambda.imag() * std::pow(site.imag(), 2) * site.real()}; */
    }

    /* template<typename func>
    std::complex<double> get_drift_term(const std::complex<double> site, func& transformer) const
    {
        std::complex<double> site_{transformer(site.real()), site.imag()};
        return mp.sigma * site_ + mp.lambda * std::pow(site_, 3) + mp.h;
    } */


    std::complex<double> get_second_order_drift_term(const std::complex<double> site)
    {
        return 3.0 * mp.lambda * std::pow(site, 2) + mp.sigma;
    }

    /* std::complex<double> get_discrete_potential(const std::complex<double> site)
    {
        return {0.5 * mp.sigma.real() * std::pow(site.real(), 2) - mp.sigma.imag() * site.imag() * site.real() +
                0.25 * mp.lambda.real() * std::pow(site.real(), 4) - 3.0/2.0 * mp.lambda.real() * std::pow(site.imag(), 2) * std::pow(site.real(), 2) +
                mp.lambda.imag() * std::pow(site.imag(), 3) * site.real() - mp.lambda.imag() * std::pow(site.real(), 3) * site.imag(),
                0.5 * mp.sigma.real() * std::pow(site.imag(), 2) + mp.sigma.imag() * site.real() * site.imag() -
                0.25 * mp.lambda.real() * std::pow(site.imag(), 4) + 3.0/2.0 * mp.lambda.real() * std::pow(site.real(), 2) * std::pow(site.imag(), 2) +
                mp.lambda.imag() * std::pow(site.real(), 3) * site.imag() - mp.lambda.imag() * std::pow(site.imag(), 3) * site.real()
        };
    } */

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

private:
    const ComplexPolynomialModelParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
