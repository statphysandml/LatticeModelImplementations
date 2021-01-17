//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_ANHARMONIC_OSCILLATOR_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_ANHARMONIC_OSCILLATOR_HPP

#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"

namespace lm_impl {
    namespace lattice_system {

        class ComplexAnharmonicOscillator;


        class ComplexAnharmonicOscillatorParameters : public LatticeModelParameters {
        public:
            explicit ComplexAnharmonicOscillatorParameters(const json params_) : LatticeModelParameters(params_),
                                                                                 dt(get_entry<double>("dt")),
                                                                                 m(get_entry<std::complex<double>>(
                                                                                         "m")),
                                                                                 omega_sq(
                                                                                         get_entry<std::complex<double>>(
                                                                                                 "omega_sq")),
                                                                                 lambda(get_entry<std::complex<double>>(
                                                                                         "lambda")) {}

            explicit ComplexAnharmonicOscillatorParameters(
                    const double dt_,
                    const std::complex<double> m_,
                    const std::complex<double> omega_sq_,
                    const std::complex<double> lambda_) :
                    ComplexAnharmonicOscillatorParameters(json{
                            {"dt",       dt_},
                            {"m",        m_},
                            {"omega_sq", omega_sq_},
                            {"lambda",   lambda_}
                    }) {}

            const static std::string name() {
                return "ComplexAnharmonicOscillator";
            }

            typedef ComplexAnharmonicOscillator Model;

        private:
            friend class ComplexAnharmonicOscillator;

            const double dt;
            const std::complex<double> m;
            const std::complex<double> omega_sq;
            const std::complex<double> lambda;
        };


        class ComplexAnharmonicOscillator : public LatticeModel<ComplexAnharmonicOscillator> {
        public:
            explicit ComplexAnharmonicOscillator(const ComplexAnharmonicOscillatorParameters &mp_) :
                    mp(mp_) {}

            std::complex<double>
            get_potential(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) {
                // neighbour[0] corresponds to the right neighbour, neighbour[1] to the left one
                return 0.5 * mp.m * (std::pow(*neighbours[0] - site, 2.0) + std::pow(site - *neighbours[1], 2.0)) /
                       mp.dt +
                       0.5 * mp.dt * mp.omega_sq * std::pow(site, 2.0) + mp.lambda / 24.0 * mp.dt * std::pow(site, 4.0);

            }

            std::complex<double>
            get_energy_per_lattice_elem(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours)
            {
                // neighbour[0] corresponds to the right neighbour, neighbour[1] to the left one
                return 0.5 * mp.m * std::pow(*neighbours[0] - site, 2.0) / mp.dt +
                       0.5 * mp.dt * mp.omega_sq * std::pow(site, 2.0) + mp.lambda / 24.0 * mp.dt * std::pow(site, 4.0);
            }

            std::complex<double>
            get_drift_term(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) {
                return mp.m * (-1.0 * (*neighbours[0] - site) + (site - *neighbours[1])) / mp.dt +
                       mp.dt * mp.omega_sq * site + mp.lambda / 6.0 * mp.dt * std::pow(site, 3.0);
            }

            std::complex<double> get_second_order_drift_term(const std::complex<double> site,
                                                             const std::vector<std::complex<double> *> neighbours) {
                return mp.m / mp.dt + mp.dt * mp.omega_sq + mp.lambda / 2.0 * mp.dt * std::pow(site, 2.0);
            }

        private:
            const ComplexAnharmonicOscillatorParameters &mp;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_ANHARMONIC_OSCILLATOR_HPP
