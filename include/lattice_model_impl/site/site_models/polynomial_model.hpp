//
// Created by lukas on 12.11.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_POLYNOMIAL_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


namespace lm_impl {
    namespace site_system {

        class PolynomialModel;


        class PolynomialModelParameters : public SiteModelParameters {
        public:
            explicit PolynomialModelParameters(const json params_) : SiteModelParameters(params_),
                                                                     sigma(get_entry<double>("sigma")),
                                                                     lambda(get_entry<double>("lambda")),
                                                                     h(get_entry<double>("h", 0.0)) {}

            explicit PolynomialModelParameters(double lambda_, double sigma_, double h_) : PolynomialModelParameters(
                    json{
                            {"lambda", lambda_},
                            {"sigma",  sigma_},
                            {"h",      h_}
                    }) {}

            static std::string name() {
                return "PolynomialModel";
            }

            typedef PolynomialModel Model;

        private:
            friend class PolynomialModel;

            const double sigma;
            const double lambda;
            const double h;
        };


        class PolynomialModel : public SiteModel<PolynomialModel> {
        public:
            explicit PolynomialModel(const PolynomialModelParameters &mp_) : mp(mp_) {}

            double get_potential(const double site) const {
                return 0.5 * mp.sigma * std::pow(site, 2) + 0.25 * mp.lambda * std::pow(site, 4) + mp.h * site;
            }

            double get_drift_term(const double site) const {
                return mp.sigma * site + mp.lambda * std::pow(site, 3) + mp.h;
            }

        private:
            const PolynomialModelParameters &mp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_POLYNOMIAL_MODEL_HPP
