//
// Created by lukas on 11.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_DUMMY_UPDATE_DYNAMICS_HPP
#define LATTICEMODELIMPLEMENTATIONS_DUMMY_UPDATE_DYNAMICS_HPP



#include "update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct DummyUpdateDynamics;

        struct DummyUpdateDynamicsParameters : UpdateDynamicsBaseParameters {
            explicit DummyUpdateDynamicsParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit DummyUpdateDynamicsParameters() : DummyUpdateDynamicsParameters(json{}) {}

            static std::string name() {
                return "DummyUpdateDynamics";
            }

            typedef DummyUpdateDynamics UpdateDynamics;
        };

        struct DummyUpdateDynamics : public UpdateDynamicsBase<DummyUpdateDynamics> {
            explicit DummyUpdateDynamics(const DummyUpdateDynamicsParameters &lp_) : lp(lp_) {}

            template<typename Site>
            void initialize_update(const Site &site) {}

            template<typename Site>
            void update(Site &site, uint measure_interval = 1) {}

            const DummyUpdateDynamicsParameters &lp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_DUMMY_UPDATE_DYNAMICS_HPP
