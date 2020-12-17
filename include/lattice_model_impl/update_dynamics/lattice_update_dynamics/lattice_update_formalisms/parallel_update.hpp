//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct ParallelUpdate;

        struct ParallelUpdateParameters : UpdateDynamicsBaseParameters {
            explicit ParallelUpdateParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit ParallelUpdateParameters() : ParallelUpdateParameters(json{}) {}

            static std::string name() {
                return "ParallelUpdate";
            }

            typedef ParallelUpdate UpdateDynamics; //  LatticeUpdate;
        };

        struct ParallelUpdate : public UpdateDynamicsBase<ParallelUpdate> {
            explicit ParallelUpdate(const ParallelUpdateParameters &lp_) : lp(lp_) {}

            template<typename Lattice>
            void initialize_update(const Lattice &lattice) {}

            template<typename Lattice>
            void update(Lattice &lattice, uint measure_interval = 1) {
                // Introduce boost?
                for (uint j = 0; j < measure_interval; j++) {
                    std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.get_size(),
                                                                             typename Lattice::SiteType(0));

                    // #pragma omp parallel for
                    for (uint i = 0; i < lattice.get_size(); i++) {
                        lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i],
                                                                  lattice.neighbours_at(i));
                    }

                    // Rewrite?
                    auto &lattice_grid = lattice.get_system_representation();
                    lattice_grid = lattice_grid_new;
                }
            }

            const ParallelUpdateParameters &lp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP
