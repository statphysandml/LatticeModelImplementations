//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_O_N_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_O_N_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"
#include "../lattice_model.hpp"

// #include "../su2.hpp"

namespace lm_impl {
    namespace lattice_system {

        class ComplexONModel;

        class ComplexONModelParameters : public LatticeModelParameters {
        public:
            explicit ComplexONModelParameters(const json params_) : LatticeModelParameters(params_),
                                                             beta(get_entry<double>("beta")),
                                                             m(get_entry<std::complex<double>>("m")),
                                                             lambda(get_entry<std::complex<double>>("lambda"))
            {}

            explicit ComplexONModelParameters(double beta_, std::complex<double> m_, std::complex<double> lambda_) : ComplexONModelParameters(json{
                    {"beta", beta_},
                    {"m", m_},
                    {"lambda", lambda_}
            }) {}

            const static std::string name() {
                return "ComplexONModel";
            }

            typedef ComplexONModel Model;

        private:
            friend class ComplexONModel;

            const double beta;
            const std::complex<double> m;
            const std::complex<double> lambda;
        };

        class ComplexONModel : public LatticeModel<ComplexONModel> {
        public:
            explicit ComplexONModel(const ComplexONModelParameters &mp_) : mp(mp_) {}

            template<typename T, typename T2=std::complex<double_t>>
            T2 get_potential(const T site, const std::vector<T*> neighbours)
            {
                double sum = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    sum += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                sum = sum - 0.5 * (2.0 * site.dim() + pow(mp.m, 2.0)) * site_sq + 0.25 * mp.lambda * pow(site_sq, 2.0);
                return  -1.0 * mp.beta * sum; // 0.5
            }

            template<typename T, typename T2=std::complex<double_t>>
            T2 get_drift_term(const T site, const std::vector<T *> neighbours) {
                return T2(0.0);
            }

        private:
            const ComplexONModelParameters &mp;
        };

        struct ComplexONModelSampler
        {
            ComplexONModelSampler(const double eps_) : eps(eps_)
            {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            T random_state() {
                T new_site(0);
                for(auto i = 0; i < new_site.dim(); i++)
                    new_site(i) += std::sqrt(2 * eps) * normal(mcmc::util::gen);
                return new_site;
            }

            template<typename T>
            T propose_state(T site) {
                return site + random_state<T>();
            }

            double get_eps() const
            {
                return eps;
            }

            const static std::string name() {
                return "ComplexONModelSampler";
            }

            const double eps;
            std::normal_distribution<double> normal;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_O_N_MODEL_HPP
