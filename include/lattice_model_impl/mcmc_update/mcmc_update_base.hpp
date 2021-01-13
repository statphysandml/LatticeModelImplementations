//
// Created by lukas on 11.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_MCMC_UPDATE_BASE_HPP
#define LATTICEMODELIMPLEMENTATIONS_MCMC_UPDATE_BASE_HPP


#include "param_helper/params.hpp"
#include "mcmc_simulation/util/random.hpp"


namespace lm_impl {
    namespace mcmc_update {

        class MCMCUpdateBaseParameters : public param_helper::params::Parameters {
        public:
            MCMCUpdateBaseParameters(const json params_) : Parameters(params_),
                                                           eps(get_entry<double>("eps", 0.0)) {}

            void write_to_file(const std::string &root_dir) {
                Parameters::write_to_file(root_dir, param_file_name());
            }

            static const std::string param_file_name() {
                return "mcmc_update_params";
            }

            const double eps;
        };

        template<typename MCMCUpdate, typename SamplerCl>
        class MCMCUpdateBase {
        public:
            MCMCUpdateBase(const double eps) : sampler(SamplerCl(eps)) {}

            template<typename Site>
            void initialize(const Site &site) {
                return mcmc_update_update().initialize_mcmc_update(site);
            }

            template<typename Site>
            void initialize_mcmc_update(const Site &site) {}

            template<typename T>
            T random_state() {
                return sampler.template random_state<T>();
            };

            template<typename T>
            T cold_state() {
                return T(0);
            };

            template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_mcmc_update_measures(const SBP &system_parameters) {
                // auto measure_names = system_parameters.get_measures();
                return std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>{};
            }


        protected:
            SamplerCl sampler;

        private:
            MCMCUpdate &mcmc_update_update() {
                return *static_cast<MCMCUpdate *>(this);
            }

            const MCMCUpdate &mcmc_update_update() const {
                return *static_cast<const MCMCUpdate *>(this);
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_MCMC_UPDATE_BASE_HPP
