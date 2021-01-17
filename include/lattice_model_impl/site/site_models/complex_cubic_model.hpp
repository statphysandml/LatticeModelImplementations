//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_COMPLEX_CUBIC_MODEL_HPP
#define MAIN_COMPLEX_CUBIC_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


namespace lm_impl {
    namespace site_system {

        class ComplexCubicModel;


        class ComplexCubicModelParameters : public SiteModelParameters {
        public:
            explicit ComplexCubicModelParameters(const json params_) : SiteModelParameters(params_) {}

            explicit ComplexCubicModelParameters() : ComplexCubicModelParameters(json{}) {}

            static std::string name() {
                return "ComplexCubicModel";
            }

            typedef ComplexCubicModel Model;

        private:
            friend class ComplexCubicModel;
        };


        class ComplexCubicModel : public SiteModel<ComplexCubicModel> {
        public:
            explicit ComplexCubicModel(const ComplexCubicModelParameters &mp_) : mp(mp_) {}

            static std::complex<double> get_drift_term(const std::complex<double> site) {
                return {-2.0 * site.real() * site.imag(), -1.0 * (std::pow(site.imag(), 2) - std::pow(site.real(), 2))};
            }

            static std::complex<double> get_second_order_drift_term(const std::complex<double> site) {
                return {-2.0 * site.imag(), 2.0 * site.real()};
            }

            static std::complex<double> get_potential(const std::complex<double> site) {
                // return std::complex<double>{0, 1} * std::pow(site, 3) / 3.0;
                return {-1.0 * std::pow(site.real(), 2) * site.imag() + std::pow(site.imag(), 3) / 3.0,
                        std::pow(site.real(), 3) / 3.0 - std::pow(site.imag(), 2) * site.real()};
            }

        private:
            const ComplexCubicModelParameters &mp;
        };

// int ComplexCubicModel::init = 0;

    }
}

#endif //MAIN_COMPLEX_CUBIC_MODEL_HPP
