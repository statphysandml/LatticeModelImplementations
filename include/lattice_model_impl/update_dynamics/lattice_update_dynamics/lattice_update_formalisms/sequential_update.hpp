//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SEQUENTIAL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SEQUENTIAL_UPDATE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct SequentialUpdate;

        struct SequentialUpdateParameters : UpdateDynamicsBaseParameters {
            explicit SequentialUpdateParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit SequentialUpdateParameters() : SequentialUpdateParameters(json{}) {}

            static std::string name() {
                return "SequentialUpdate";
            }

            typedef SequentialUpdate UpdateDynamics; //  LatticeUpdate;
        };

        struct SequentialUpdate : public UpdateDynamicsBase<SequentialUpdate> {
            explicit SequentialUpdate(const SequentialUpdateParameters &lp_) : lp(lp_) {}

            template<typename Lattice>
            void initialize_update(const Lattice &lattice) {
                uniint = std::uniform_int_distribution<int>(0, lattice.get_size() - 1);
            }

            template<typename Lattice>
            void update(Lattice &lattice, uint measure_interval = 1) {
                for (size_t k = 0; k < measure_interval; k++) {
                    for (uint j = 0; j < lattice.get_size(); j++) {
                        int i = uniint(mcmc::util::gen);
                        lattice[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i],
                                                         lattice.neighbours_at(i));
                        // const double K = std::fabs(update_formalism->estimate_drift_term(lattice[i], lattice.neighbours_at[i]));
                    }
                }
            }

            const SequentialUpdateParameters &lp;
            std::uniform_int_distribution<int> uniint;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SEQUENTIAL_UPDATE_HPP
