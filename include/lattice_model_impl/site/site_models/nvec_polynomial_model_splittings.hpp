//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMICAL_MODEL_SPLITTINGS_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMICAL_MODEL_SPLITTINGS_HPP

#include "../../site_representations/representations/nvec.hpp"

namespace site_models {
    namespace nvec_polynomial_model {
        namespace splittings {
            struct NVec2a
            {
                NVec2a(const double sigma, const double lambda, const double h) : sigma(sigma), lambda(lambda), h(h)
                {}

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    return {0.25 * lambda * std::pow(site.reduce(), 4), 0.5 * sigma * std::pow(site.reduce(), 2) + h * site.reduce()};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    return {lambda * std::pow(site.orig(), 3.0) + 3.0 * lambda * site.orig() * std::pow(site(1), 2.0) + 1.0 * sigma * site.orig(),
                            lambda * std::pow(site(1), 3.0) + 3.0 * lambda * std::pow(site.orig(), 2.0) * site(1) + 1.0 * sigma * site(1)};

                }

                typedef NVec<double, 2>::Ttype Ttype;

                const double sigma, lambda, h;
            };

            struct NVec2b
            {
                NVec2b(const double sigma, const double lambda, const double h) : sigma(sigma), lambda(lambda), h(h)
                {}

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    return {0.5 * sigma * std::pow(site.reduce(), 2), h * site.reduce(), 0.25 * lambda * std::pow(site.reduce(), 4)};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    return {2.0 * lambda * std::pow(site.reduce(), 3) + h, -1.0 * sigma * site.reduce() + h, 2.0* sigma * site.reduce() -1.0 * lambda * std::pow(site.reduce(), 3)};

                }

                typedef NVec<double, 2>::Ttype Ttype;

                const double sigma, lambda, h;
            };

            struct NVec6a
            {
                NVec6a(const double sigma, const double lambda, const double h) : sigma(sigma), lambda(lambda), h(h)
                {}

                std::vector<double> get_potential(const NVec<double, 6> &site) const
                {
                    return {0.25 * lambda * std::pow(site.reduce(), 4) + 0.5 * sigma * std::pow(site.reduce(), 2) + h * site.reduce(),
                            0.0, 0.0, 0.0, 0.0, 0.0};

                }

                std::vector<double> get_drift_term(const NVec<double, 6> &site) const
                {
                    return {
                            lambda * std::pow(site.orig(), 3.0),
                            lambda * std::pow(site(1), 3.0),
                            -1.0 * sigma * site(1) + 2.0 * sigma * site.orig(),
                            3.0 * lambda * site.orig() * std::pow(site(1), 2.0),
                            3.0 * lambda * std::pow(site.orig(), 2.0) * site(1),
                            2.0 * sigma * site(1) - 1.0 * sigma * site.orig()
                    };

                }

                typedef NVec<double, 6>::Ttype Ttype;

                const double sigma, lambda, h;
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMICAL_MODEL_SPLITTINGS_HPP
