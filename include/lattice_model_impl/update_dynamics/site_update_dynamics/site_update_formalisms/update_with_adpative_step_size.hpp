//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct UpdateWithAdpativeStepsize;

        struct UpdateWithAdpativeStepsizeParameters : UpdateDynamicsBaseParameters {
            explicit UpdateWithAdpativeStepsizeParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {
                thermalization_steps = get_entry<int>("thermalization_steps", 2000);
            }

            static std::string name() {
                return "UpdateWithAdpativeStepsize";
            }

            typedef UpdateWithAdpativeStepsize UpdateDynamics; // SiteUpdate;

            int thermalization_steps;
        };

        struct UpdateWithAdpativeStepsize : public UpdateDynamicsBase<UpdateWithAdpativeStepsize> {
            explicit UpdateWithAdpativeStepsize(const UpdateWithAdpativeStepsizeParameters &sp_) : sp(sp_) {}

            template<typename Site>
            void initialize_update(const Site &site) {
                thermalization_counter = 0;
            }

            template<typename Site>
            void update(Site &lattice, uint measure_interval = 1) {
                // (this->*function_ptr_on_update_func) (lattice);
                if (thermalization_counter > sp.thermalization_steps)
                    parallel_update_with_adpative_stepsize(lattice, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(lattice, measure_interval);
            }

            template<typename Site>
            void thermalization_phase_with_adpative_stepsize(Site &site, uint measure_interval = 1) {
                // ToDo Adapt measure interval
                if (measure_interval > 1)
                    std::cout << "Thermalization phase with adaptive stepsize needs to be completed" << std::endl;

                const double KMax = std::fabs(
                        site.get_update_formalism().estimate_drift_term(site.get_system_representation()));

                if (thermalization_counter < sp.thermalization_steps) {
                    KExpectation += KMax;
                    thermalization_counter++;
                } else if (thermalization_counter == sp.thermalization_steps) {
                    KExpectation /= sp.thermalization_steps;
                    thermalization_counter = 0;
                }

                site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                       site.get_system_representation());
            }


            template<typename Site>
            void parallel_update_with_adpative_stepsize(Site &site, uint measure_interval = 1) {
                for (auto k = 0; k < measure_interval; k++) {
                    const double KMax = std::fabs(
                            site.get_update_formalism()->estimate_drift_term(site.get_system_representation()));
                    site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                           site.get_system_representation(), KMax,
                                                                           KExpectation);
                }
            }

            const UpdateWithAdpativeStepsizeParameters &sp;

            // For adaptive step size
            double KExpectation = 0;
            int thermalization_counter = 0;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
