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

            explicit UpdateWithAdpativeStepsizeParameters(int thermalization_steps_) :
                    UpdateWithAdpativeStepsizeParameters(json {
                        {"thermalization_steps", thermalization_steps_}})
            {}

            static std::string name() {
                return "UpdateWithAdpativeStepsize";
            }

            typedef UpdateWithAdpativeStepsize UpdateDynamics;

            int thermalization_steps;
        };

        struct UpdateWithAdpativeStepsize : public UpdateDynamicsBase<UpdateWithAdpativeStepsize> {
            explicit UpdateWithAdpativeStepsize(const UpdateWithAdpativeStepsizeParameters &sp_) : sp(sp_) {}

            template<typename Site>
            void initialize_update(const Site &site) {
                thermalized = false;
            }

            template<typename Site>
            void update(Site &lattice, uint measure_interval = 1) {
                if (thermalized)
                    parallel_update_with_adpative_stepsize(lattice, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(lattice);
            }

            template<typename Site>
            void thermalization_phase_with_adpative_stepsize(Site &site) {
                for (int k = 0; k < sp.thermalization_steps; k++) {
                    KExpectation += std::fabs(site.get_update_formalism().estimate_drift_term(site.get_system_representation()));
                    site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                           site.get_system_representation());
                }

                KExpectation /= sp.thermalization_steps;
                thermalized = true;
            }


            template<typename Site>
            void parallel_update_with_adpative_stepsize(Site &site, uint measure_interval = 1) {
                for (uint k = 0; k < measure_interval; k++) {
                    const double KMax = std::fabs(site.get_update_formalism().estimate_drift_term(site.get_system_representation()));
                    site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                           site.get_system_representation(), KMax,
                                                                           KExpectation);
                }
            }

            const UpdateWithAdpativeStepsizeParameters &sp;

            // For adaptive step size
            double KExpectation = 0;
            bool thermalized = 0;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
