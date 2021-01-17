//
// Created by lukas on 12.01.21.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_LINK_LATTICE_MODEL_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_LINK_LATTICE_MODEL_MEASURES_HPP

#include "mcmc_simulation/measure_policy.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/params.hpp"
#include "param_helper/json.hpp"


namespace lm_impl {
    namespace util {
        namespace link_lattice_system_model_measures {
            template<typename SB>
            struct MeasurePolyakovLoopPolicy: public mcmc::common_measures::MeasurePolicy< SB > {
            public:
                MeasurePolyakovLoopPolicy(const std::vector<int> dimensions_, const uint elem_per_site_) : dimensions(dimensions_), elem_per_site(elem_per_site_)
                {}

                std::string measure(const SB &system) override {
                    // Hasn't been tested!!
                    std::complex<double> polyakov_loop = 0;
                    for(auto i = 0; i < system.size(); i += elem_per_site * dimensions[0])
                    {
                        typename SB::SiteType group_elem = system[i];
                        for(auto tau = elem_per_site; tau < elem_per_site * dimensions[0]; tau += elem_per_site)
                            group_elem = group_elem * system[i + tau];
                        polyakov_loop += group_elem.trace();
                    }
                    return std::to_string(polyakov_loop);
                }

                std::string name()
                {
                    return "PolyakovLoop";
                }

                const std::vector<int> dimensions; // Different dimensions
                const uint elem_per_site; // Number of elements per site
            };

            template<typename SB>
            struct MeasureAveragePlaquetteActionPolicy : public mcmc::common_measures::MeasurePolicy<SB> {
            public:
                MeasureAveragePlaquetteActionPolicy(const int dimension_) : n_plaquettes_per_link((dimension_ - 1) * 2)
                {}

                // Needs to be devided by the inverse temperature in a postprocessing step
                std::string measure(const SB &system) override {
                    return std::to_string(system.energy() / n_plaquettes_per_link);
                }

                std::string name() {
                    return "AveragePlaquetteAction";
                }

                const int n_plaquettes_per_link;
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_LINK_LATTICE_MODEL_MEASURES_HPP
