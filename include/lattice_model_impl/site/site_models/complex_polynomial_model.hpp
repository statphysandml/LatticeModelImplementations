//
// Created by lukas on 02.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


namespace lm_impl {
    namespace site_system {

        class ComplexPolynomialModel;


        class ComplexPolynomialModelParameters : public SiteModelParameters {
        public:
            explicit ComplexPolynomialModelParameters(const json params_) : SiteModelParameters(params_),
                                                                            sigma(std::complex<double> {get_entry<double>("sigma_real", 0.0),
                                                                                    get_entry<double>("sigma_imag", 0.0)}),
                                                                            lambda(std::complex<double> {get_entry<double>("lambda_real", 0.0),
                                                                                    get_entry<double>("lambda_imag", 0.0)}),
                                                                            h(std::complex<double> {get_entry<double>("h_real", 0.0),
                                                                                    get_entry<double>("h_imag", 0.0)})
            {}

            explicit ComplexPolynomialModelParameters(double lambda_real_, double lambda_imag_, double sigma_real_,
                                                      double sigma_imag_, double h_real_, double h_imag_)
                    : ComplexPolynomialModelParameters(json{
                    {"lambda_real", lambda_real_},
                    {"lambda_imag", lambda_imag_},
                    {"sigma_real", sigma_real_},
                    {"sigma_imag", sigma_imag_},
                    {"h_real", h_real_},
                    {"h_imag", h_imag_},
            }) {}

            static std::string name() {
                return "ComplexPolynomialModel";
            }

            typedef ComplexPolynomialModel Model;

        private:
            friend class ComplexPolynomialModel;

            const std::complex<double> sigma;
            const std::complex<double> lambda;
            const std::complex<double> h;
        };


        class ComplexPolynomialModel : public SiteModel<ComplexPolynomialModel> {
        public:
            explicit ComplexPolynomialModel(const ComplexPolynomialModelParameters &mp_) : mp(mp_) {}

            std::complex<double> get_potential(const std::complex<double> site) const {
                return 0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site;
            }

            std::complex<double> get_drift_term(const std::complex<double> site) const {
                return mp.sigma * site + mp.lambda * std::pow(site, 3) + mp.h;
            }

            std::complex<double> get_second_order_drift_term(const std::complex<double> site) const {
                return 3.0 * mp.lambda * std::pow(site, 2) + mp.sigma;
            }

        private:
            const ComplexPolynomialModelParameters &mp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
