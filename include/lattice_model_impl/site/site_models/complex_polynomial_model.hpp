//
// Created by lukas on 02.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"



class ComplexPolynomialModel;


class ComplexPolynomialModelParameters : public SiteModelParameters {
public:
    explicit ComplexPolynomialModelParameters(const json params_) : SiteModelParameters(params_),
                                                                               sigma_real(get_value_by_key< double >("sigma_real", 0.0)),
                                                                               sigma_imag(get_value_by_key< double >("sigma_imag", 0.0)),
                                                                               sigma(get_value_by_key< std::complex<double> >("sigma", {sigma_real, sigma_imag})),
                                                                               lambda(get_value_by_key< std::complex<double> >("lambda")),
                                                                               h(get_value_by_key< std::complex<double> >("h", 0.0))
    {}

    explicit ComplexPolynomialModelParameters(std::complex<double> lambda_, std::complex<double> sigma_, std::complex<double> h_) : ComplexPolynomialModelParameters(json {
            {"lambda", lambda_},
            {"sigma", sigma_},
            {"h", h_}
    })
    {}

    explicit ComplexPolynomialModelParameters(std::complex<double> lambda_, double sigma_real_, double sigma_imag_, std::complex<double> h_) : ComplexPolynomialModelParameters(json {
            {"lambda", lambda_},
            {"sigma_real", sigma_real_},
            {"sigma_imag", sigma_imag_},
            {"sigma", {sigma_real_, sigma_imag_}},
            {"h", h_}
    })
    {}

    static std::string name() {
        return "ComplexPolynomialModel";
    }

    typedef ComplexPolynomialModel Model;

private:
    friend class ComplexPolynomialModel;

    const double sigma_real;
    const double sigma_imag;
    const std::complex<double> sigma;
    const std::complex<double> lambda;
    const std::complex<double> h;
};


class ComplexPolynomialModel : public SiteModel< ComplexPolynomialModel >
{
public:
    explicit ComplexPolynomialModel(const ComplexPolynomialModelParameters &mp_) : mp(mp_) {}

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

    std::complex<double> get_second_order_drift_term(const std::complex<double> site) const
    {
        return 3.0 * mp.lambda * std::pow(site, 2) + mp.sigma;
    }

    // [
    // For Cobrid monte Carlo Algorithms

    double get_potential(double site) const
    {
        return (0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site).real();
    }

    double get_imag_potential(double site) const
    {
        return (0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site).imag();
    }

    double get_drift_term(const double site) const
    {
        return (mp.sigma * site + mp.lambda * std::pow(site, 3) + mp.h).real();
    }

    double get_imag_drift_term(const double site) const
    {
        return (mp.sigma * site + mp.lambda * std::pow(site, 3) + mp.h).imag();
        /* return  {mp.sigma.real() * site.real() - mp.sigma.imag() * site.imag() +
                    mp.lambda.real() * std::pow(site.real(), 3) - 3.0 * mp.lambda.real() * std::pow(site.imag(), 2) * site.real() +
                    mp.lambda.imag() * std::pow(site.imag(), 3) - 3.0 * mp.lambda.imag() * std::pow(site.real(), 2) * site.imag(),
                 mp.sigma.imag() * site.real() + mp.sigma.real() * site.imag() -
                    mp.lambda.real() * std::pow(site.imag(), 3) + 3.0 * mp.lambda.real() * std::pow(site.real(), 2) * site.imag() +
                    mp.lambda.imag() * std::pow(site.real(), 3) - 3.0 * mp.lambda.imag() * std::pow(site.imag(), 2) * site.real()}; */
    }

    // ]

    // [
    // For Complex Monte Carlo

    double get_real_const_imag_drift_term(const std::complex<double> site)
    {
        /* auto derivative_real = mp.h + std::pow(site, 3) * mp.lambda + site * mp.sigma;
        auto derivative_imag = std::complex<double>{0.0, 1.0} * (mp.h + std::pow(site, 3) * mp.lambda + site * mp.sigma); */
        return get_drift_term(site).imag();
    }

    double get_imag_const_imag_drift_term(const std::complex<double> site)
    {
        /* auto derivative_real = mp.h + std::pow(site, 3) * mp.lambda + site * mp.sigma;
        auto derivative_imag = std::complex<double>{0.0, 1.0} * (mp.h + std::pow(site, 3) * mp.lambda + site * mp.sigma); */
        return get_drift_term(site).real();
    }

    // ]

    // [
    // For Complex Hybrid Monte Carlo

    std::complex<double> get_complex_hybrid_drift_term(const std::complex<double> site)
    {
        return mp.lambda.imag() * std::pow(site, 3.0) + mp.sigma.imag() * site;
    }

    std::complex<double> get_complex_hybrid_potential(const std::complex<double> site) const
    {
        return 0.5 * mp.sigma.imag() * std::pow(site, 2) + 0.25 * mp.lambda.imag() * std::pow(site, 4) + mp.h.imag() * site;
    }

    std::complex<double> get_real_complex_hybrid_potential(const std::complex<double> site) const
    {
        return 0.5 * mp.sigma.real() * std::pow(site, 2) + 0.25 * mp.lambda.real() * std::pow(site, 4) + mp.h.real() * site;
    }

    // ]

    /* template<typename func>
    std::complex<double> get_potential(const std::complex<double> site, func& transformer) const
    {
        std::complex<double> site_{transformer(site.real()), site.imag()};
        return 0.5 * mp.sigma * std::pow(site_, 2) + 0.25 * mp.lambda * std::pow(site_, 4) + mp.h * site_;
    } */

    /* template<typename func>
    std::complex<double> get_drift_term(const std::complex<double> site, func& transformer) const
    {
        std::complex<double> site_{transformer(site.real()), site.imag()};
        return mp.sigma * site_ + mp.lambda * std::pow(site_, 3) + mp.h;
    } */

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

private:
    const ComplexPolynomialModelParameters &mp;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
