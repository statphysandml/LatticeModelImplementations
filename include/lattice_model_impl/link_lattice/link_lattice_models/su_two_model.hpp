//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_SU_TWO_MODEL_HPP
#define MAIN_SU_TWO_MODEL_HPP


#include "../link_lattice_model.hpp"


namespace lm_impl {
    namespace link_lattice_system {


        class SUTwoModel;

        class SUTwoModelParameters : public LinkLatticeModelParameters {
        public:
            explicit SUTwoModelParameters(const json params_) : LinkLatticeModelParameters(params_),
                                                                beta(get_entry<double>("beta"))
                                                                {}

            explicit SUTwoModelParameters(double beta_) : SUTwoModelParameters(json{
                    {"beta", beta_}
            }) {}

            const static std::string name() {
                return "SUTwoModel";
            }

            static uint N() {
                return 2;
            };

            typedef SUTwoModel Model;

        private:
            friend class SUTwoModel;

            const double beta;
        };


        class SUTwoModel : public LinkLatticeModel<SUTwoModel> {
        public:
            explicit SUTwoModel(const SUTwoModelParameters &mp_) : mp(mp_) {}

            // Corresponds to the contribution of the link to the action: S[U] = \beta/N tr[neighbours.size() / 3 * identity - U * A]
            // See Gattringer - Quantum Chromodynamics on the Lattice Section 4.1.4
            template<typename T, typename T2=double_t>
            T2 get_potential(const T link, const std::vector<T *> neighbours) {
                T A = calc_A(neighbours);
                return mp.beta / SUTwoModelParameters::N() * (neighbours.size() / 3 * SUTwoModelParameters::N() -
                                                              std::real((link * A).trace()));
            }

            template<typename T, typename T2=double_t>
            T2 get_drift_term(const T link, const std::vector<T *> neighbours) {
                return T2(0.0);
            }

        private:
            const SUTwoModelParameters &mp;
        };

        struct SUTwoModelSampler
        {
            SUTwoModelSampler(const double eps_) : eps(eps_)
            {}

            template<typename T>
            T random_state() {
                return T("random");
            }

            template<typename T>
            T propose_state(T link) {
                return link * T(eps);
            }

            double get_eps() const
            {
                return eps;
            }

            const static std::string name() {
                return "SUTwoModelSampler";
            }

            const double eps;
        };
    }
}

#endif //MAIN_SU_TWO_MODEL_HPP
