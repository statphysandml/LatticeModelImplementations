//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct ParallelUpdateWithAdpativeStepsize;

        struct ParallelUpdateWithAdpativeStepsizeParameters : UpdateDynamicsBaseParameters {
            explicit ParallelUpdateWithAdpativeStepsizeParameters(const json params_) : UpdateDynamicsBaseParameters(
                    params_) {
                thermalization_steps = get_entry<int>("thermalization_steps", 2000);
            }

            explicit ParallelUpdateWithAdpativeStepsizeParameters(int thermalization_steps_) :
                    ParallelUpdateWithAdpativeStepsizeParameters(json {
                            {"thermalization_steps", thermalization_steps_}})
            {}

            static std::string name() {
                return "ParallelUpdateWithAdpativeStepsize";
            }

            typedef ParallelUpdateWithAdpativeStepsize UpdateDynamics;

            int thermalization_steps;
        };

        struct ParallelUpdateWithAdpativeStepsize : public UpdateDynamicsBase<ParallelUpdateWithAdpativeStepsize> {
            explicit ParallelUpdateWithAdpativeStepsize(const ParallelUpdateWithAdpativeStepsizeParameters &lp_) : lp(
                    lp_) {}

            template<typename Lattice>
            void initialize_update(const Lattice &lattice) {
                thermalized = false;
            }

            template<typename Lattice>
            void update(Lattice &lattice, uint measure_interval = 1) {
                if (thermalized)
                    parallel_update_with_adpative_stepsize(lattice, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(lattice);
            }

            template<typename Lattice>
            void thermalization_phase_with_adpative_stepsize(Lattice &lattice, uint measure_interval = 1) {
                // static_assert(detail::is_updateable<T, typename UpdateFormalismParameters::MCMCUpdate>::value, "is not estimate_drift_term");

                for (int k = 0; k < lp.thermalization_steps; k++) {

                    double KMax = 0;
                    for (uint i = 0; i < lattice.get_size(); i++) {
                        const double K = std::fabs(
                                lattice.get_update_formalism().estimate_drift_term(lattice[i],
                                                                                   lattice.neighbours_at(i)));
                        if (K > KMax)
                            KMax = K;
                    }
                    KExpectation += KMax;

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

                KExpectation /= lp.thermalization_steps;
                thermalized = true;
            }


            template<typename Lattice>
            void parallel_update_with_adpative_stepsize(Lattice &lattice, uint measure_interval = 1) {
                for (uint j = 0; j < measure_interval; j++) {

                    double KMax = 0;
                    for (uint i = 0; i < lattice.get_size(); i++) {
                        const double K = std::fabs(
                                lattice.get_update_formalism().estimate_drift_term(lattice[i],
                                                                                   lattice.neighbours_at(i)));
                        if (K > KMax)
                            KMax = K;
                    }

                    std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.get_size(),
                                                                             typename Lattice::SiteType(0));

                    // #pragma omp parallel for
                    for (uint i = 0; i < lattice.get_size(); i++) {
                        lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i],
                                                                  lattice.neighbours_at(i), KMax, KExpectation);
                    }

                    // Rewrite?
                    auto &lattice_grid = lattice.get_system_representation();
                    lattice_grid = lattice_grid_new;
                }
            }

            const ParallelUpdateWithAdpativeStepsizeParameters &lp;

            // For adaptive step size
            double KExpectation = 0;
            bool thermalized = 0;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
