//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_MEMORY_SIMPLE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_MEMORY_SIMPLE_UPDATE_HPP


#include "../../../update_dynamics/update_dynamics_base.hpp"
#include "../../../util/measures/complex_monte_carlo_measures.hpp"


namespace lm_impl {
    namespace update_dynamics {

        template<typename SiteType>
        struct MemorySiteSimpleUpdate;

        template<typename SiteType>
        struct MemorySiteSimpleUpdateParameters : UpdateDynamicsBaseParameters {
            explicit MemorySiteSimpleUpdateParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {}

            explicit MemorySiteSimpleUpdateParameters() : MemorySiteSimpleUpdateParameters(json{}) {}

            static std::string name() {
                return "MemorySiteSimpleUpdate";
            }

            typedef MemorySiteSimpleUpdate<SiteType> UpdateDynamics;
        };

        template<typename SiteType>
        struct MemorySiteSimpleUpdate : public UpdateDynamicsBase<MemorySiteSimpleUpdate<SiteType>> {
            explicit MemorySiteSimpleUpdate(const MemorySiteSimpleUpdateParameters<SiteType> &lp_) : lp(lp_) {}

            template<typename Site>
            void initialize_update(const Site &site) {
                site_system_params = site.get_params_json();
            }

            template<typename Site>
            void update(Site &site, uint measure_interval = 1) {
                for (size_t k = 0; k < measure_interval - 1; k++) {
                    site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                           site.get_system_representation());
                }
                previous_site = site.get_system_representation();
                site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                       site.get_system_representation());
            }

            template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_update_dynamics_measures(const SBP &system_parameters) {
                auto measure_names = system_parameters.get_measures();

                std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>> measures{};
                for (auto &measure_name :  measure_names)
                    if (measure_name == "DetailedBalanceAccuracy")
                        measures.push_back(
                                std::make_unique<util::complex_monte_carlo_measures::MeasureDetailedBalanceAccuracyPolicy<SB>>(
                                        site_system_params, previous_site));
                    else if (measure_name == "AbsoluteDetailedBalanceAccuracy")
                        measures.push_back(
                                std::make_unique<util::complex_monte_carlo_measures::MeasureAbsoluteDetailedBalanceAccuracyPolicy<SB>>(
                                        site_system_params, previous_site));
                    else if (measure_name == "RealStepSize")
                        measures.push_back(
                                std::make_unique<util::complex_monte_carlo_measures::MeasureRealStepSizePolicy<SB>>(
                                        previous_site));
                    else if (measure_name == "ImagStepSize")
                        measures.push_back(
                                std::make_unique<util::complex_monte_carlo_measures::MeasureImagStepSizePolicy<SB>>(
                                        previous_site));
                    else if (measure_name == "ComplexStepSize")
                        measures.push_back(
                                std::make_unique<util::complex_monte_carlo_measures::MeasureComplexStepSizePolicy<SB>>(
                                        previous_site));
                return measures;
            }

            const MemorySiteSimpleUpdateParameters<SiteType> &lp;

            SiteType previous_site;
            json site_system_params;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_MEMORY_SIMPLE_UPDATE_HPP
