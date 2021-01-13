//
// Created by lukas on 09.10.19.
//

#ifndef MAIN_MODEL_HPP
#define MAIN_MODEL_HPP


#include "param_helper/params.hpp"
#include "mcmc_simulation/measure_policy.hpp"


namespace lm_impl {
    namespace lattice_system {

        class LatticeModelParameters : public param_helper::params::Parameters {
        public:
            using Parameters::Parameters;

            const void write_to_file(const std::string &root_dir) {
                Parameters::write_to_file(root_dir, param_file_name());
            }

            static const std::string param_file_name() {
                return "model_params";
            }
        };


        template<typename ModelCL>
        class LatticeModel {
        public:
            template<typename T>
            T normalize(T state) {
                return state;
            }

            template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_model_measures(const SBP &system_parameters) {
                // auto measure_names = system_parameters.get_measures();
                return std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>{};
            }

        private:
            ModelCL &lattice_model() {
                return *static_cast<ModelCL *>(this);
            }
        };

    }
}

#endif //MAIN_MODEL_HPP
