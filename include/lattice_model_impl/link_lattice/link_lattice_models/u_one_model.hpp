//
// Created by lukas on 04.11.19.
//

#ifndef MAIN_U_ONE_MODEL_HPP
#define MAIN_U_ONE_MODEL_HPP

#include "../link_lattice_model.hpp"

namespace lm_impl {
    namespace link_lattice_system {


        class UOneModel;

        class UOneModelParameters : public LinkLatticeModelParameters {
        public:
            explicit UOneModelParameters(const json params_) : LinkLatticeModelParameters(params_),
                                                               beta(get_entry<double>("beta")),
                                                               eps(get_entry<double>("eps"))
            //mu(complex_to_json(get_entry("mu")))
            {}

            explicit UOneModelParameters(double beta_, double eps_) : UOneModelParameters(json{
                    {"beta", beta_},
                    {"eps",  eps_}
            }) {}

            const static std::string name() {
                return "UOneModel";
            }

            static uint N() {
                return 1;
            };

            typedef UOneModel Model;

        private:
            friend class UOneModel;

            const double beta;
            const double eps;
        };


        class UOneModel : public LinkLatticeModel<UOneModel> {
        public:
            explicit UOneModel(const UOneModelParameters &mp_) : mp(mp_) {}

            template<typename T>
            T random_state() {
                return T("random");
            }

            template<typename T>
            T propose_state(T site) {
                return site * T(mp.eps);
            }

            template<typename T, typename T2=double_t>
            T2 get_potential(const T site, const std::vector<T *> neighbours) {
                T A = calc_A(neighbours);
                // (lat.dim()-1)*2-std::real(1.0/N*(lat(n*lat.dim()+mu)*calc_A(lat, n, mu, both_orientations)).trace())
                return mp.beta / mp.N() * (neighbours.size() / 3.0 -
                                           std::real((site * A).trace())); // Additional 1/N missing in second term??
            }

        private:
            const UOneModelParameters &mp;
        };

    }
}

#endif //MAIN_U_ONE_MODEL_HPP
