//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMICAL_MODEL_SPLITTINGS_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMICAL_MODEL_SPLITTINGS_HPP

#include "mcmc_simulation/sampler/gaussian_sampler.hpp"

namespace site_models {
    namespace nvec_polynomial_model {
        namespace sampler {
            struct GaussianNVecSampler : GaussianSampler
            {
                using GaussianSampler::GaussianSampler;

                template<typename T>
                T random_state()
                {
                    return T(GaussianSampler::random_state<typename T::Ttype>());
                }

                template<typename T>
                T propose_state(T site)
                {
                    return site + T(GaussianSampler::random_state<typename T::Ttype>());
                }

                const static std::string name() {
                    return "GaussianNVecSampler";
                }
            };
        }

        namespace splittings {
            // Inherit!!
            struct NVecSplitting
            {
                NVecSplitting(const double sigma, const double lambda, const double h) : sigma(sigma), lambda(lambda), h(h)
                {}

                const double sigma, lambda, h;
            };

            struct NVec2a : NVecSplitting
            {
                using NVecSplitting::NVecSplitting;

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    // 1/4*(x^4 + 6 x^2*y^2 + 4 x*y^3) + 1/2*y^2, 1/4*(4*x^3*y + y^4) + 1/2*(x^2 + 2x*y)
                    return {0.25 * lambda * (std::pow(site(0), 4.0) +
                                             6.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             4.0 * site(0) * std::pow(site(1), 3.0)) +
                            0.5 * sigma * std::pow(site(1), 2.0) +
                            h * site(1),
                            0.25 * lambda * (4.0 * std::pow(site(0), 3.0) * site(1) +
                                             std::pow(site(1), 4.0)) +
                            0.5 * sigma * (std::pow(site(0), 2.0) +
                                           2.0 * site(0) * site(1)) +
                            h * site(0)};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    // x^3 + 3 x * y^2 + y^3, 3 * x^2 * y + x + y
                    return {lambda * (std::pow(site(0), 3.0) +
                                      3.0 * site(0) * std::pow(site(1), 2.0) +
                                      std::pow(site(1), 3.0)),
                            lambda * 3.0 * std::pow(site(0), 2.0) * site(1) +
                            sigma * (site(0) + site(1)) +
                            h};
                }

                static std::string name() {
                    return "NVec2a";
                }

                static std::string get_drift_term_in_python_code()
                {
                    return "drift_term = lambda x, y: (np.power(x, 3.0) + 3.0 * x * np.power(y, 2.0) + np.power(y, 3.0), 3.0 * np.power(x, 2.0) * y + x + y)";
                }

                typedef NVec<double, 2>::Ttype Ttype;
            };

            struct NVec2b : NVecSplitting
            {
                using NVecSplitting::NVecSplitting;

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    // 1/4*(-x^4 + 6 x^2*y^2 + 4 x*y^3) + 1/2*y^2, 1/4*(2*x^4 + 4*x^3*y + y^4) + 1/2*(x^2 + 2 x*y)
                    return {0.25 * lambda * (-1.0 * std::pow(site(0), 4.0) +
                                             6.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             4.0 * site(0) * std::pow(site(1), 3.0)) +
                            0.5 * sigma * std::pow(site(1), 2.0) +
                            h * site(1),
                            0.25 * lambda * (2.0 * std::pow(site(0), 4.0) +
                                             4.0 * std::pow(site(0), 3.0) * site(1) +
                                             std::pow(site(1), 4.0)) +
                            0.5 * sigma * (std::pow(site(0), 2.0) +
                                     2.0 * site(0) * site(1)) +
                            h * site(0)};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    // -x^3 + 3 x * y^2 + y^3, 2* x^3 + 3 * x^2 * y + x + y
                    return {lambda * (-1.0 * std::pow(site(0), 3.0) +
                                      3.0 * site(0) * std::pow(site(1), 2.0) +
                                      std::pow(site(1), 3.0)),
                            lambda * (2.0 * std::pow(site(0), 3.0) +
                                      3.0 * std::pow(site(0), 2.0) * site(1)) +
                            sigma * (site(0) + site(1)) +
                            h};
                }

                static std::string name() {
                    return "NVec2b";
                }

                static std::string get_drift_term_in_python_code()
                {
                    return "drift_term = lambda x, y: (-1.0 * np.power(x, 3.0) + 3.0 * x * np.power(y, 2.0) + np.power(y, 3.0), 2.0 * np.power(x, 3.0) + 3.0 * np.power(x, 2.0) * y + x + y)";
                }

                typedef NVec<double, 2>::Ttype Ttype;
            };

            struct NVec2c : NVecSplitting
            {
                using NVecSplitting::NVecSplitting;

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    // 1/4*(-x^4 - 4*x^3*y + 4 x^2*y^2 + 4 x*y^3) + 1/2*y^2, 1/4*(2*x^4 + 8*x^3*y + 2 x^2*y^2 + y^4) + 1/2*(x^2 + 2 x*y)
                    return {0.25 * lambda * (-1.0 * std::pow(site(0), 4.0) -
                                             -4.0 * std::pow(site(0), 3.0) * site(1) +
                                             4.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             4.0 * site(0) * std::pow(site(1), 3.0)) +
                            0.5 * sigma * std::pow(site(1), 2.0) +
                            h * site(1),
                            0.25 * lambda * (2.0 * std::pow(site(0), 4.0) +
                                             8.0 * std::pow(site(0), 3.0) * site(1) +
                                             2.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             std::pow(site(1), 4.0)) +
                            0.5 * sigma * (std::pow(site(0), 2.0) +
                                           2.0 * site(0) * site(1)) +
                            h * site(0)};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    // -x^3 + 2 x * y^2 - 3 * x^2 * y + y^3, 2* x^3 + x * y^2 + 6 * x^2 * y + x + y
                    return {lambda * (-1.0 * std::pow(site(0), 3.0) +
                                      2.0 * site(0) * std::pow(site(1), 2.0) -
                                      3.0 * std::pow(site(0), 2.0) * site(1) +
                                      std::pow(site(1), 3.0)),
                            lambda * (2.0 * std::pow(site(0), 3.0) +
                                      1.0 * site(0) * std::pow(site(1), 2.0) +
                                      6.0 * std::pow(site(0), 2.0) * site(1)) +
                            sigma * (site(0) + site(1)) +
                            h};
                }

                static std::string name() {
                    return "NVec2c";
                }

                static std::string get_drift_term_in_python_code()
                {
                    return "drift_term = lambda x, y: (-1.0 * np.power(x, 3.0) + 2.0 * x * np.power(y, 2.0) - 3.0 * np.power(x, 2.0) * y + np.power(y, 3.0), 2.0 * np.power(x, 3.0) + x * np.power(y, 2.0) + 6.0 * np.power(x, 2.0) * y + x + y)";
                }

                typedef NVec<double, 2>::Ttype Ttype;
            };

            struct NVec2d : NVecSplitting
            {
                using NVecSplitting::NVecSplitting;

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    // 1/4*(-x^4 + 4 x^3*y - 2 x^2*y^2 + 4 x*y^3) + 1/2*y^2, 1/4*(2 x^4 + 8 x^2*y^2 + y^4) + 1/2*(x^2 + 2 x*y)
                    return {0.25 * lambda * (-1.0 * std::pow(site(0), 4.0) -
                                             4.0 * std::pow(site(0), 3.0) * site(1) -
                                             2.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             4.0 * site(0) * std::pow(site(1), 3.0)) +
                            0.5 * sigma * std::pow(site(1), 2.0) +
                            h * site(1),
                            0.25 * lambda * (2.0 * std::pow(site(0), 4.0) +
                                             8.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             std::pow(site(1), 4.0)) +
                            0.5 * sigma * (std::pow(site(0), 2.0) +
                                           2.0 * site(0) * site(1)) +
                            h * site(0)};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    // -x^3 + 3 x^2 * y - x * y^2+ y^3, 2* x^3 + 4 x * y^2 + x + y
                    return {lambda * (-1.0 * std::pow(site(0), 3.0) +
                                      3.0 * std::pow(site(0), 2.0) * site(1) -
                                      1.0 * site(0) * std::pow(site(1), 2.0) +
                                      std::pow(site(1), 3.0)),
                            lambda * (2.0 * std::pow(site(0), 3.0) +
                                      4.0 * site(0) * std::pow(site(1), 2.0)) +
                            sigma * (site(0) + site(1)) +
                            h};
                }

                static std::string name() {
                    return "NVec2d";
                }

                static std::string get_drift_term_in_python_code()
                {
                    return "drift_term = lambda x, y: (-1.0 * np.power(x, 3.0) + 3.0 * np.power(x, 2.0) * y - 1.0 * x * np.power(y, 2.0) + np.power(y, 3.0), 2.0 * np.power(x, 3.0) + 4.0 * x * np.power(y, 2.0) + x + y)";
                }

                typedef NVec<double, 2>::Ttype Ttype;
            };

            struct NVec2x
            {
                NVec2x(const double sigma, const double lambda, const double h) : sigma(sigma), lambda(lambda), h(h)
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

                static std::string name() {
                    return "NVec2x";
                }

                typedef NVec<double, 2>::Ttype Ttype;

                const double sigma, lambda, h;
            };

            struct NVec3
            {
                NVec3(const double sigma, const double lambda, const double h) : sigma(sigma), lambda(lambda), h(h)
                {}

                std::vector<double> get_potential(const NVec<double, 3> &site) const
                {
                    return {0.5 * sigma * std::pow(site.reduce(), 2), h * site.reduce(), 0.25 * lambda * std::pow(site.reduce(), 4)};

                }

                std::vector<double> get_drift_term(const NVec<double, 3> &site) const
                {
                    return {2.0 * lambda * std::pow(site.reduce(), 3) + h, -1.0 * sigma * site.reduce() + h, 2.0* sigma * site.reduce() -1.0 * lambda * std::pow(site.reduce(), 3)};

                }

                static std::string name() {
                    return "NVec3";
                }

                typedef NVec<double, 3>::Ttype Ttype;

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

                static std::string name() {
                    return "NVec6a";
                }

                typedef NVec<double, 6>::Ttype Ttype;

                const double sigma, lambda, h;
            };

            struct NVec2VariableSplitting : NVecSplitting
            {
                NVec2VariableSplitting(const double sigma, const double lambda, const double h) : NVecSplitting(sigma, lambda, h)
                {
                    ft = std::vector<double> (7, 0.0);
                    normal = std::normal_distribution<double> (0,0.01);
                    uniint = std::uniform_int_distribution<int>(0, 6);
                }

                std::vector<double> ft;
                std::normal_distribution<double> normal;
                std::uniform_int_distribution<int> uniint;

                void update_ft()
                {
                    auto ft_index = uniint(gen);
                    ft[ft_index] += normal(gen);
                }

                std::vector<double> get_potential(const NVec<double, 2> &site) const
                {
                    // 1/4*(-x^4 + 4 x^3*y - 2 x^2*y^2 + 4 x*y^3) + 1/2*y^2, 1/4*(2 x^4 + 8 x^2*y^2 + y^4) + 1/2*(x^2 + 2 x*y)
                    return {0.25 * lambda * (std::pow(site(0), 4.0) -
                                             4.0 * std::pow(site(0), 3.0) * site(1) -
                                             2.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             4.0 * site(0) * std::pow(site(1), 3.0)) +
                            0.5 * sigma * std::pow(site(1), 2.0) +
                            h * site(1),
                            0.25 * lambda * (2.0 * std::pow(site(0), 4.0) +
                                             8.0 * std::pow(site(0), 2.0) * std::pow(site(1), 2.0) +
                                             std::pow(site(1), 4.0)) +
                            0.5 * sigma * (std::pow(site(0), 2.0) +
                                           2.0 * site(0) * site(1)) +
                            h * site(0)};

                }

                std::vector<double> get_drift_term(const NVec<double, 2> &site) const
                {
                    // -x^3 + 3 x^2 * y - x * y^2+ y^3, 2* x^3 + 4 x * y^2 + x + y
                    return {lambda * ((1.0 - ft[0]) * std::pow(site(0), 3.0) +
                                      (3.0 - ft[1]) * std::pow(site(0), 2.0) * site(1) +
                                      (3.0 - ft[2]) * site(0) * std::pow(site(1), 2.0) +
                                      (1.0 - ft[3]) * std::pow(site(1), 3.0)) +
                            sigma * ((1.0 - ft[4]) * site(0) + (1.0 - ft[5]) * site(1)) +
                            h * (1.0 - ft[6]),
                            lambda * (ft[0] * std::pow(site(0), 3.0) +
                                      ft[1] * std::pow(site(0), 2.0) * site(1) +
                                      ft[2] * site(0) * std::pow(site(1), 2.0) +
                                      ft[3] * std::pow(site(1), 3.0)) +
                            sigma * (ft[4] * site(0) + ft[5] * site(1)) +
                            h * ft[6]};
                }

                static std::string name() {
                    return "NVec2VariableSplitting";
                }

                static std::string get_drift_term_in_python_code()
                {
                    return "drift_term = lambda x, y: (-1.0 * np.power(x, 3.0) + 3.0 * np.power(x, 2.0) * y - 1.0 * x * np.power(y, 2.0) + np.power(y, 3.0), 2.0 * np.power(x, 3.0) + 4.0 * x * np.power(y, 2.0) + x + y)";
                }

                typedef NVec<double, 2>::Ttype Ttype;
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMICAL_MODEL_SPLITTINGS_HPP
