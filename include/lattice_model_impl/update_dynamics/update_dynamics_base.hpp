//
// Created by lukas on 11.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_UPDATE_DYNAMICS_BASE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UPDATE_DYNAMICS_BASE_HPP


#include "param_helper/params.hpp"
#include "mcmc_simulation/measure_policy.hpp"

namespace lm_impl {
    namespace update_dynamics {

// ToDo: auto update_lattice_site(UpdateFormalism& f, Args&... args) with the &, the neighbours are destroyed after an update -> still true?

// https://www.grimm-jaud.de/index.php/blog/perfect-forwarding
// https://stackoverflow.com/questions/3582001/what-are-the-main-purposes-of-using-stdforward-and-which-problems-it-solves
        template<typename UpdateFormalism, typename... Args>
        auto update_lattice_site(UpdateFormalism &f, Args &&... args) {
            return f(std::forward<Args>(args)...);
        }

        template<typename UpdateFormalism, typename... Args>
        void global_lattice_update(UpdateFormalism &f, Args &&... args) {
            f(std::forward<Args>(args)...);
        }


        class UpdateDynamicsBaseParameters : public param_helper::params::Parameters {
        public:
            using Parameters::Parameters;

            const void write_to_file(const std::string &root_dir) {
                Parameters::write_to_file(root_dir, param_file_name());
            }

            static const std::string param_file_name() {
                return "update_dynamics_params";
            }
        };


        template<typename Derived>
        class UpdateDynamicsBase {
        public:
            template<typename System>
            void initialize(const System &system) {
                return system_update().initialize_update(system);
            }

            template<typename System>
            void operator()(System &system, uint measure_interval = 1) {
                return system_update().update(system, measure_interval);
            }

            template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_update_dynamics_measures(const SBP &system_parameters) {
                // auto measure_names = system_parameters.get_measures();
                return std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>{};
            }

        private:
            Derived &system_update() {
                return *static_cast<Derived *>(this);
            }

            const Derived &system_update() const {
                return *static_cast<const Derived *>(this);
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_UPDATE_DYNAMICS_BASE_HPP
