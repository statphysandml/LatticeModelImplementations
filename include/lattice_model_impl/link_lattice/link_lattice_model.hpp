//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_LINK_LATTICE_MODEL_HPP
#define MAIN_LINK_LATTICE_MODEL_HPP

#include "../lattice/lattice_model.hpp"

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
            template<typename SB>
        struct MeasurePolyakovLoopPolicy: public mcmc::common_measures::MeasurePolicy< SB > {
            public:
                MeasurePolyakovLoopPolicy(const std::vector<int> &dimensions_, const uint &elem_per_site_) : dimensions(dimensions_), elem_per_site(elem_per_site_)
                {}

                std::string measure(SB &system) override {
                    std::complex<double> polyakov_loop = 0;
                    for(auto i = 0; i < system.size(); i += elem_per_site * dimensions[0])
                    {
                        typename SB::SiteType group_elem = system[i];
                        for(auto tau = elem_per_site; tau < elem_per_site * dimensions[0]; tau += elem_per_site)
                            group_elem = group_elem * system[i + tau];
                        polyakov_loop += group_elem.trace();
                    }
                    //return std::to_string(TypePolicy<double>::imagv(system.energy()));
                    return std::to_string(polyakov_loop.real()) + " " + std::to_string(polyakov_loop.imag());
                }

                std::string name()
                {
                    return "PolyakovLoop";
                }

                const std::vector<int> &dimensions; // Different dimensions
                const uint &elem_per_site; // Number of elements per site
            };

            template<typename SB, typename SBP>
            mcmc::common_measures::MeasurePolicy<SB>* model_measure_factory(const std::string& measure, const SBP& system_base_parameters) {
                if(measure == "PolyakovLoop")
                    return new MeasurePolyakovLoopPolicy<SB>(system_base_parameters.get_dimensions(), system_base_parameters.get_elems_per_site());
                else
                    return nullptr;
            }
        protected:
            template<typename T>
            T calc_A (const std::vector<T*> neighbours, bool both_orientations=true) {
                T A("null");
                for(auto i = 0; i < neighbours.size(); i += 6)
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
