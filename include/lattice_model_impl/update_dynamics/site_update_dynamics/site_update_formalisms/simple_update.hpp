//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct SiteSimpleUpdate;

        struct SiteSimpleUpdateParameters : UpdateDynamicsBaseParameters {
            explicit SiteSimpleUpdateParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit SiteSimpleUpdateParameters() : SiteSimpleUpdateParameters(json{}) {}

            static std::string name() {
                return "SiteSimpleUpdate";
            }

            typedef SiteSimpleUpdate UpdateDynamics; // SiteUpdate;
        };

        struct SiteSimpleUpdate : public UpdateDynamicsBase<SiteSimpleUpdate> {
            explicit SiteSimpleUpdate(const SiteSimpleUpdateParameters &lp_) : lp(lp_) {}

            template<typename Site>
            void initialize_update(const Site &site) {}

            template<typename Site>
            void update(Site &site, uint measure_interval = 1) {
                for (uint k = 0; k < measure_interval; k++) {
                    site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                           site.get_system_representation());
                }
            }

            const SiteSimpleUpdateParameters &lp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP
