//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct GlobalLatticeUpdate;

        struct GlobalLatticeUpdateParameters : UpdateDynamicsBaseParameters {
            explicit GlobalLatticeUpdateParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit GlobalLatticeUpdateParameters() : GlobalLatticeUpdateParameters(json{}) {}

            static std::string name() {
                return "GlobalLatticeUpdate";
            }

            typedef GlobalLatticeUpdate UpdateDynamics; // LatticeUpdate;
        };

        struct GlobalLatticeUpdate : public UpdateDynamicsBase<GlobalLatticeUpdate> {
            explicit GlobalLatticeUpdate(const GlobalLatticeUpdateParameters &lp_) : lp(lp_) {}

            template<typename Lattice>
            void initialize_update(const Lattice &lattice) {}

            template<typename Lattice>
            void update(Lattice &lattice, uint measure_interval = 1) {
                for (uint k = 0; k < measure_interval; k++) {
                    global_lattice_update(lattice.get_update_formalism(), lattice);
                }
            }

            const GlobalLatticeUpdateParameters &lp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP
