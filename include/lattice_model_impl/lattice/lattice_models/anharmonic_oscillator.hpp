//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_ANHARMONIC_OSCILLATOR_HPP
#define LATTICEMODELIMPLEMENTATIONS_ANHARMONIC_OSCILLATOR_HPP

#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"

template<typename SamplerCl>
class AnharmonicOscillator;

template<typename SamplerCl>
class AnharmonicOscillatorParameters : public LatticeModelParameters {
public:
    explicit AnharmonicOscillatorParameters(const json params_) : LatticeModelParameters(params_),
                                                                  dt(get_value_by_key<double>("dt")),
                                                                  m(get_value_by_key<double>("m")),
                                                                  omega_sq(get_value_by_key<double>("omega_sq")),
                                                                  lambda(get_value_by_key<double>("lambda")),
                                                                  eps(get_value_by_key<double>("eps", 0.1))
    {}

    explicit AnharmonicOscillatorParameters(
            const double eps_,
            const double dt_,
            const double m_,
            const double omega_sq_,
            const double lambda_) :
            AnharmonicOscillatorParameters(json{
                    {"eps", eps_},
                    {"dt", dt_},
                    {"m", m_},
                    {"omega_sq", omega_sq_},
                    {"lambda", lambda_}
            })
    {}

    const static std::string name() {
        return "AnharmonicOscillator";
    }

    typedef AnharmonicOscillator<SamplerCl> Model;

private:
    friend class AnharmonicOscillator<SamplerCl>;

    const double eps;
    const double dt;
    const double m;
    const double omega_sq;
    const double lambda;
};

template<typename SamplerCl>
class AnharmonicOscillator : public LatticeModel< AnharmonicOscillator<SamplerCl> >
{
public:
    explicit AnharmonicOscillator(const AnharmonicOscillatorParameters<SamplerCl> &mp_) :
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

    template<typename T>
    T get_drift_term(const T site, const std::vector<T*> neighbours)
    {
        return mp.m * (-1.0 * (*neighbours[0] - site) + (site - *neighbours[1])) / mp.dt +
               mp.dt * mp.omega_sq * site + mp.lambda / 6.0 * mp.dt * std::pow(site, 3.0);
    }

    template<typename T>
    T get_second_order_drift_term(const T site, const std::vector<T*> neighbours)
    {
        return mp.m / mp.dt + mp.dt * mp.omega_sq + mp.lambda / 2.0 * mp.dt * std::pow(site, 2.0);
    }

    template<typename T>
    T get_potential(const T site, const std::vector<T*> neighbours)
    {
        // neighbour[0] corresponds to the right neighbour, neighbour[1] to the left one
        return 0.5 * mp.m * (std::pow(*neighbours[0] - site, 2.0) + std::pow(site - *neighbours[1], 2.0)) / mp.dt +
               0.5 * mp.dt * mp.omega_sq * std::pow(site, 2.0) + mp.lambda / 24.0 * mp.dt * std::pow(site, 4.0);

    }

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

private:
    const AnharmonicOscillatorParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

#endif //LATTICEMODELIMPLEMENTATIONS_ANHARMONIC_OSCILLATOR_HPP
