//
// Created by lukas on 31.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP


#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"


namespace lm_impl {
    namespace lattice_system {

        template<typename SB>
        struct MeasureXYMagnetization : public mcmc::common_measures::MeasurePolicy<SB> {
        public:
            std::string measure(const SB &system) override {
                auto sum_cos = decltype(system[0]){0};
                auto sum_sin = decltype(system[0]){0};

                for (auto i = 0; i < system.size(); i++) {
                    sum_cos += std::cos(system[i]);
                    sum_sin += std::sin(system[i]);
                }

                return std::to_string(std::pow(sum_cos / system.size(), 2) + std::pow(sum_sin / system.size(), 2));
            }

            std::string name() {
                return "XYMagnetization";
            }
        };


        class XYModel;


        class XYModelParameters : public LatticeModelParameters {
        public:
            explicit XYModelParameters(const json params_) : LatticeModelParameters(params_),
                                                             beta(get_entry<double>("beta", 0.5)),
                                                             J(get_entry<double>("J", 1.0)),
                                                             h(get_entry<double>("h", 0.0)) {}

            explicit XYModelParameters(double beta_, double J_, double h_, double eps_) : XYModelParameters(json{
                    {"beta", beta_},
                    {"J",    J_},
                    {"h",    h_}
            }) {}

            const static std::string name() {
                return "XYModel";
            }

            typedef XYModel Model;

        private:
            friend class XYModel;

            const double beta;
            const double J;
            const double h;
        };


        class XYModel : public LatticeModel<XYModel> {
        public:
            explicit XYModel(const XYModelParameters &mp_) : mp(mp_) {}

            template<typename T>
            T normalize(T state) {
                state = std::fmod(state, 2 * M_PI);
                if (state < 0) {
                    state += 2 * M_PI;
                }
                return state;
            }

            template<typename T>
            T get_potential(const T site, const std::vector<T *> neighbours) {
                double S = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S += mp.J * std::cos(site - *neighbours[i]) + mp.h * std::cos(site);
                }
                return -1.0 * mp.beta * S;
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours)
            {
                double S = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S += mp.J * std::cos(site - *neighbours[i]) + mp.h * std::cos(site);
                }
                return -1.0 * mp.beta * S;
            }

            template<typename T>
            T get_drift_term(const T site, const std::vector<T *> neighbours) {
                double S = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S += mp.J * std::sin(site - *neighbours[i]) +
                         mp.J * std::sin(site - *neighbours[i + 1]) + mp.h * std::cos(site);
                }
                return mp.beta * S;
            }

            template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_model_measures(const SBP &system_parameters) {
                auto measure_names = system_parameters.get_measures();

                std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>> measures{};
                for (auto &measure_name :  measure_names)
                    if (measure_name == "XYMagnetization")
                        measures.push_back(std::make_unique<MeasureXYMagnetization <SB>>());
                return measures;
            }

        private:
            const XYModelParameters &mp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP
