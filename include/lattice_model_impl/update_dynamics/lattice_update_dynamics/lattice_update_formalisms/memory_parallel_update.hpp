//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_MEMORY_PARALLEL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_MEMORY_PARALLEL_UPDATE_HPP

#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        template<typename SiteType>
        struct MemoryParallelUpdate;

        template<typename SiteType>
        struct MemoryParallelUpdateParameters : UpdateDynamicsBaseParameters {
            explicit MemoryParallelUpdateParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit MemoryParallelUpdateParameters() : MemoryParallelUpdateParameters(json{}) {}

            static std::string name() {
                return "MemoryParallelUpdate";
            }

            typedef MemoryParallelUpdate<SiteType> UpdateDynamics;
        };

        template<typename SiteType>
        struct MemoryParallelUpdate : public UpdateDynamicsBase<MemoryParallelUpdate<SiteType>> {
            explicit MemoryParallelUpdate(const MemoryParallelUpdateParameters<SiteType> &lp_) : lp(lp_) {}

            template<typename Lattice>
            void initialize_update(const Lattice &lattice) {}

            template<typename Lattice>
            void update(Lattice &lattice, uint measure_interval = 1) {
                // ToDo: Introduce boost!
                for (auto j = 0; j < measure_interval; j++) {
                    std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.get_size(),
                                                                             typename Lattice::SiteType(0));

                    // #pragma omp parallel for
                    for (uint i = 0; i < lattice.get_size(); i++) {
                        lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i],
                                                                  lattice.neighbours_at(i));
                    }

                    if (j == measure_interval - 1)
                        previous_lattice = lattice.get_syste_representation();

                    // ToDo: Rewrite?
                    auto &lattice_grid = lattice.get_system_representation();
                    lattice_grid = lattice_grid_new;
                }
            }


            // ToDo: Can be defined only once outside of the class

            /* template<typename SB>
            struct MeasureDetailedBalanceAccuracyPolicy : public mcmc::common_measures::MeasurePolicy<SB> {
            public:
                explicit MeasureDetailedBalanceAccuracyPolicy(std::vector<SiteType> &prev_lattice_) : prev_lattice(
                        prev_lattice_) {}

                std::string measure(const SB &system) override {
                    auto a = prev_lattice;
                    return std::to_string(system.energy());
                }

                std::string name() {
                    return "DetailedBalanceAccuracy";
                }

                std::vector<SiteType> &prev_lattice;
            };

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_update_dynamics_measures(const json &measure_names) {
                std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>> measures{};
                for (auto &measure_name :  measure_names)
                    if (measure_name == "DetailedBalanceAccuracy")
                        measures.push_back(
                                std::make_unique<MeasureDetailedBalanceAccuracyPolicy<SB>>(previous_lattice));
                return measures;
            } */

            const MemoryParallelUpdateParameters<SiteType> &lp;
            std::vector<SiteType> previous_lattice;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_MEMORY_PARALLEL_UPDATE_HPP
