//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP


#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"


namespace lm_impl {
    namespace lattice_system {

        class ComplexXYModel;


        class ComplexXYModelParameters : public LatticeModelParameters {
        public:
            explicit ComplexXYModelParameters(const json params_) : LatticeModelParameters(params_),
                                                                    beta(get_entry<double>("beta")),
                                                                    mu(get_entry<double>("mu"))
            {}

            explicit ComplexXYModelParameters(double beta_, double mu_) : ComplexXYModelParameters(
                    json{
                            {"beta", beta_},
                            {"mu",   mu_}
                    }) {}

            const static std::string name() {
                return "ComplexXYModel";
            }

            typedef ComplexXYModel Model;

        private:
            friend class ComplexXYModel;

            const double beta;
            const double mu;
        };


        class ComplexXYModel : public LatticeModel<ComplexXYModel> {
        public:
            explicit ComplexXYModel(const ComplexXYModelParameters &mp_) :
                    mp(mp_) {
            }

            std::complex<double> normalize(std::complex<double> state) {
                state.real(normalize(state.real()));
                return state;
            }

            double normalize(double state) {
                state = std::fmod(state, 2 * M_PI);
                if (state < 0) {
                    state += 2 * M_PI;
                }
                return state;
            }

            std::complex<double>
            get_potential(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) {
                double S_re = 0;
                double S_im = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S_re += std::cos(site.real() - neighbours[i]->real()) *
                            std::cosh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                            std::cos(neighbours[i + 1]->real() - site.real()) *
                            std::cosh(neighbours[i + 1]->imag() - site.imag() - mp.mu * int(i == 0));
                    S_im += std::sin(site.real() - neighbours[i]->real()) *
                            std::sinh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                            std::sin(neighbours[i + 1]->real() - site.real()) *
                            std::sinh(neighbours[i + 1]->imag() - site.imag() - mp.mu * int(i == 0));
                }
                return {-1.0 * mp.beta * S_re, mp.beta * S_im};
            }

            std::complex<double>
            get_energy_per_lattice_elem(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) {
                double S_re = 0;
                double S_im = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S_re += std::cos(site.real() - neighbours[i]->real()) *
                            std::cosh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0));
                    S_im += std::sin(site.real() - neighbours[i]->real()) *
                            std::sinh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0));
                }
                return {-1.0 * mp.beta * S_re, mp.beta * S_im};
            }

            std::complex<double>
            get_drift_term(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) {
                double S_re = 0;
                double S_im = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S_re += std::sin(site.real() - neighbours[i]->real()) *
                            std::cosh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                            std::sin(site.real() - neighbours[i + 1]->real()) *
                            std::cosh(site.imag() - neighbours[i + 1]->imag() + mp.mu * int(i == 0));
                    S_im += std::cos(site.real() - neighbours[i]->real()) *
                            std::sinh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                            std::cos(site.real() - neighbours[i + 1]->real()) *
                            std::sinh(site.imag() - neighbours[i + 1]->imag() + mp.mu * int(i == 0));
                }
                return mp.beta * std::complex<double>(S_re, S_im);
            }

        private:
            const ComplexXYModelParameters &mp;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP
