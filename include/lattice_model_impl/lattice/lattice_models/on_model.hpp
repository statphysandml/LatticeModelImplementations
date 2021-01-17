//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP

#include "mcmc_simulation/util/random.hpp"
#include "../lattice_model.hpp"

namespace lm_impl {
    namespace lattice_system {

        class ONModel;

        class ONModelParameters : public LatticeModelParameters {
        public:
            explicit ONModelParameters(const json params_) : LatticeModelParameters(params_),
                                                                beta(get_entry<double>("beta")),
                                                                kappa(get_entry<double>("kappa")),
                                                                lambda(get_entry<double>("lambda"))
            {}

            explicit ONModelParameters(double beta_, double m_, double lambda_) : ONModelParameters(json{
                    {"beta", beta_},
                    {"kappa", m_},
                    {"lambda", lambda_}
            }) {}

            const static std::string name() {
                return "ONModel";
            }

            typedef ONModel Model;

        private:
            friend class ONModel;

            const double beta;
            const double kappa;
            const double lambda;
        };

        class ONModel : public LatticeModel<ONModel> {
        public:
            explicit ONModel(const ONModelParameters &mp_) : mp(mp_) {}

            // According to equation (3.16) from Smit QFT Lattice
            template<typename T, typename T2=double_t>
            T2 get_potential(const T site, const std::vector<T*> neighbours)
            {
                double potential = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    potential += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                potential = 2.0 * mp.kappa * potential - site_sq - mp.lambda * pow(site_sq - 1.0, 2.0);
                return -1.0 * mp.beta * potential;
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours)
            {
                double potential = 0;
                for(size_t i = 0; i < neighbours.size(); i += 2) {
                    potential += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                potential = 2.0 * mp.kappa * potential - site_sq - mp.lambda * pow(site_sq - 1.0, 2.0);
                return -1.0 * mp.beta * potential;
            }

            template<typename T, typename T2=double_t>
            T2 get_drift_term(const T site, const std::vector<T *> neighbours) {
                return T2(0.0);
            }

        private:
            const ONModelParameters &mp;
        };

        // Overload GaussianSampler instead?

        struct ONModelSampler
        {
            ONModelSampler(const double eps_) : eps(eps_)
            {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            T random_state() {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++)
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
                return "ONModelSampler";
            }

            const double eps;
            std::normal_distribution<double> normal;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP
