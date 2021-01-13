//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_LINK_LATTICE_MODEL_HPP
#define MAIN_LINK_LATTICE_MODEL_HPP

#include "../lattice/lattice_model.hpp"

#include "../util/measures/link_lattice_model_measures.hpp"

namespace lm_impl {
    namespace link_lattice_system {

        class LinkLatticeModelParameters : public lattice_system::LatticeModelParameters {
        public:
            using LatticeModelParameters::LatticeModelParameters;
        };


        template <typename ModelCL>
        class LinkLatticeModel : public lattice_system::LatticeModel<ModelCL>
        {
        public:
            template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_model_measures(const SBP &system_parameters) {
                auto measure_names = system_parameters.get_measures();
                std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>> link_lattice_measures{};
                for (auto &measure_name :  measure_names)
                    if (measure_name == "PolyakovLoop")
                        link_lattice_measures.push_back(std::make_unique<util::link_lattice_system_model_measures::MeasurePolyakovLoopPolicy<SB>>(
                                system_parameters.get_dimensions(), system_parameters.get_elems_per_site()));
                    else if (measure_name == "AveragePlaquetteAction")
                        link_lattice_measures.push_back(std::make_unique<util::link_lattice_system_model_measures::MeasureAveragePlaquetteActionPolicy<SB>>(
                                system_parameters.get_dimension()));
                return link_lattice_measures;
            }

        protected:
            template<typename T>
            T calc_A (const std::vector<T*> neighbours, bool both_orientations=true) {
                T A("null");
                for(uint i = 0; i < neighbours.size(); i += 6)
                {
                    A += (*neighbours[i]) * T(*neighbours[i + 1]).adjungate() * T(*neighbours[i + 2]).adjungate();
                    if(both_orientations)
                        A += T(*neighbours[i + 3]).adjungate() * T(*neighbours[ i + 4]).adjungate() * (*neighbours[i + 5]);
                }
                return A;
            }
        };

    }
}


#endif //MAIN_LINK_LATTICE_MODEL_HPP
