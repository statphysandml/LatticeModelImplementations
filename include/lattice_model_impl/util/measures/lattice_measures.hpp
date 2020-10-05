//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP

#include "mcmc_simulation/measure_policy.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/params.hpp"
#include "param_helper/json.hpp"


namespace lattice_model_measures {
    template<typename SB>
    struct MeasureEnergyPolicy : public common_measures::MeasurePolicy<SB> {
    public:
        std::string measure(const SB &system) override {
            return std::to_string(common_measures::TypePolicy<typename SB::SiteType>::realv(system.energy()));
        }

        std::string name() {
            return "Energy";
        }
    };

    template<typename SB>
    struct MeasureEnergyImagPolicy : public common_measures::MeasurePolicy<SB> {
    public:
        std::string measure(const SB &system) override {
            return std::to_string(common_measures::TypePolicy<typename SB::SiteType>::imagv(system.energy()));
        }

        std::string name() {
            return "EnergyImag";
        }
    };

    template<typename SB>
    struct MeasureDriftTermPolicy : public common_measures::MeasurePolicy<SB> {
    public:
        std::string measure(const SB &system) override {
            auto drift_term = system.drift_term();
            return std::to_string(common_measures::TypePolicy<decltype(drift_term)>::realv(drift_term)) + " " +
                   std::to_string(common_measures::TypePolicy<decltype(drift_term)>::imagv(drift_term));
        }

        std::string name() {
            return "Drift";
        }
    };

    template<typename SB>
    struct MeasureWilsonActionPolicy : public common_measures::MeasurePolicy<SB> {
    public:
        std::string measure(const SB &system) override {
            return std::to_string(common_measures::TypePolicy<typename SB::SiteType>::realv(system.energy()));
        }

        std::string name() {
            return "WilsonAction";
        }
    };
}

#endif //LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP
