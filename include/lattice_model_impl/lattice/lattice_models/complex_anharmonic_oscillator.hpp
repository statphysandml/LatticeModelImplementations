//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_ANHARMONIC_OSCILLATOR_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_ANHARMONIC_OSCILLATOR_HPP

#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"

template<typename SamplerCl>
class ComplexAnharmonicOscillator;

template<typename SamplerCl>
class ComplexAnharmonicOscillatorParameters : public LatticeModelParameters {
public:
    explicit ComplexAnharmonicOscillatorParameters(const json params_) : LatticeModelParameters(params_),
                                                                         dt(get_value_by_key<double>("dt")),
                                                                         m(get_value_by_key<std::complex<double>>("m")),
                                                                         omega_sq(get_value_by_key<std::complex<double>>("omega_sq")),
                                                                         lambda(get_value_by_key<std::complex<double>>("lambda")),
                                                                         eps(get_value_by_key<double>("eps", 0.1))
    {}

    explicit ComplexAnharmonicOscillatorParameters(
            const double dt_,
            const std::complex<double> m_,
            const std::complex<double> omega_sq_,
            const std::complex<double> lambda_,
            const double eps_) :
            ComplexAnharmonicOscillatorParameters(json{
                    {"eps", eps_},
                    {"dt", dt_},
                    {"m", m_},
                    {"omega_sq", omega_sq_},
                    {"lambda", lambda_}
            })
    {}

    const static std::string name() {
        return "ComplexAnharmonicOscillator";
    }

    typedef ComplexAnharmonicOscillator<SamplerCl> Model;

private:
    friend class ComplexAnharmonicOscillator<SamplerCl>;

    const double eps;
    const double dt;
    const std::complex<double> m;
    const std::complex<double> omega_sq;
    const std::complex<double> lambda;
};

template<typename SamplerCl>
class ComplexAnharmonicOscillator : public LatticeModel< ComplexAnharmonicOscillator<SamplerCl> >
{
public:
    explicit ComplexAnharmonicOscillator(const ComplexAnharmonicOscillatorParameters<SamplerCl> &mp_) :
            mp(mp_), sampler(SamplerCl(mp.eps))
    {}

    template<typename T>
    T random_state()
    {
        return sampler.template random_state<T>();
    }

    template<typename T>
    T propose_state(T site)
    {
        return sampler.template propose_state<T>(site);
    }

    std::complex<double> get_potential(const std::complex<double> site, const std::vector<std::complex<double>*> neighbours)
    {
        // neighbour[0] corresponds to the right neighbour, neighbour[1] to the left one
        return 0.5 * mp.m * (std::pow(*neighbours[0] - site, 2.0) + std::pow(site - *neighbours[1], 2.0)) / mp.dt +
               0.5 * mp.dt * mp.omega_sq * std::pow(site, 2.0) + mp.lambda / 24.0 * mp.dt * std::pow(site, 4.0);

    }

    std::complex<double> get_drift_term(const std::complex<double> site, const std::vector<std::complex<double>*> neighbours)
    {
        return mp.m * (-1.0 * (*neighbours[0] - site) + (site - *neighbours[1])) / mp.dt +
               mp.dt * mp.omega_sq * site + mp.lambda / 6.0 * mp.dt * std::pow(site, 3.0);
    }

    std::complex<double> get_second_order_drift_term(const std::complex<double> site, const std::vector<std::complex<double>*> neighbours)
    {
        return mp.m / mp.dt + mp.dt * mp.omega_sq + mp.lambda / 2.0 * mp.dt * std::pow(site, 2.0);
    }

    // [
    // For Cobrid Monte Carlo algorithms

    double get_potential(const double site, const std::vector<double*> neighbours)
    {
        return (0.5 * mp.m * (std::pow(*neighbours[0] - site, 2.0) + std::pow(site - *neighbours[1], 2.0)) / mp.dt +
               0.5 * mp.dt * mp.omega_sq * std::pow(site, 2.0) + mp.lambda / 24.0 * mp.dt * std::pow(site, 4.0)).real();
    }

    double get_imag_potential(const double site, const std::vector<double*> neighbours)
    {
        return (0.5 * mp.m * (std::pow(*neighbours[0] - site, 2.0) + std::pow(site - *neighbours[1], 2.0)) / mp.dt +
                0.5 * mp.dt * mp.omega_sq * std::pow(site, 2.0) + mp.lambda / 24.0 * mp.dt * std::pow(site, 4.0)).imag();
    }

    double get_drift_term(const double site, const std::vector<double*> neighbours)
    {
        return (mp.m * (-1.0 * (*neighbours[0] - site) + (site - *neighbours[1])) / mp.dt +
                mp.dt * mp.omega_sq * site + mp.lambda / 6.0 * mp.dt * std::pow(site, 3.0)).real();
    }

    double get_imag_drift_term(const double site, const std::vector<double*> neighbours)
    {
        return (mp.m * (-1.0 * (*neighbours[0] - site) + (site - *neighbours[1])) / mp.dt +
                mp.dt * mp.omega_sq * site + mp.lambda / 6.0 * mp.dt * std::pow(site, 3.0)).imag();
    }

    // ]

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

private:
    const ComplexAnharmonicOscillatorParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_ANHARMONIC_OSCILLATOR_HPP
