//
// Created by lukas on 12.01.21.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_ON_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_ON_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"
#include "../lattice_model.hpp"

// #include "../su2.hpp"

namespace lm_impl {
    namespace lattice_system {

        class ComplexONModel;

        class ComplexONModelParameters : public LatticeModelParameters {
        public:
            explicit ComplexONModelParameters(const json params_) :
                LatticeModelParameters(params_),
                beta(get_entry<double>("beta")),
                kappa(std::complex<double> {get_entry<double>("kappa_real", 0.0),
                        get_entry<double>("kappa_imag", 0.0)}),
                lambda(std::complex<double> {get_entry<double>("lambda_real", 0.0),
                        get_entry<double>("lambda_imag", 0.0)})
            {}

            explicit ComplexONModelParameters(double beta_, double kappa_real_, double kappa_imag_,
                                              double lambda_real_, double lambda_imag_) : ComplexONModelParameters(json{
                    {"beta", beta_},
                    {"kappa_real", kappa_real_},
                    {"kappa_imag", kappa_imag_},
                    {"lambda_real", lambda_real_},
                    {"lambda_imag", lambda_imag_}
            }) {}

            const static std::string name() {
                return "ComplexONModel";
            }

            typedef ComplexONModel Model;

        private:
            friend class ComplexONModel;

            const double beta;
            const std::complex<double> kappa;
            const std::complex<double> lambda;
        };

        class ComplexONModel : public LatticeModel<ComplexONModel> {
        public:
            explicit ComplexONModel(const ComplexONModelParameters &mp_) : mp(mp_) {}

            // According to equation (3.16) from Smit QFT Lattice
            template<typename T, typename T2=std::complex<double_t>>
            T2 get_potential(const T site, const std::vector<T*> neighbours)
            {
                std::complex<double> potential = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    potential += site * (*neighbours[i]);
                }
                std::complex<double> site_sq = site * site;
                potential = 2.0 * mp.kappa * potential - site_sq - mp.lambda * pow(site_sq - 1.0, 2.0);
                return -1.0 * mp.beta * potential;
            }

            template<typename T>
            T get_drift_term(const T site, const std::vector<T *> neighbours) {
                T drift_term(0);
                for(size_t i = 0; i < neighbours.size(); i++) {
                    drift_term += (*neighbours[i]);
                }
                drift_term = 2.0 * mp.kappa * drift_term + (-2.0) * site + (-4.0) * site * mp.lambda * (site * site - 1.0);
                return -1.0 * mp.beta * drift_term;
            }

        private:
            const ComplexONModelParameters &mp;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_ON_MODEL_HPP
