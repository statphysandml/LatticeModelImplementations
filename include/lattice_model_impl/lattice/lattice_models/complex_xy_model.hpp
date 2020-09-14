//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP


#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"


template<typename SamplerCl>
class ComplexXYModel;


template<typename SamplerCl>
class ComplexXYModelParameters : public LatticeModelParameters {
public:
    explicit ComplexXYModelParameters(const json params_) : LatticeModelParameters(params_),
                                                     beta(get_value_by_key<double>("beta")),
                                                     mu(get_value_by_key<double>("mu")),
                                                     eps(get_value_by_key<double>("eps", 0.1))
    //mu(complex_to_json(get_value_by_key("mu")))
    {}

    explicit ComplexXYModelParameters(double beta_, double mu_, double eps_=0.1) : ComplexXYModelParameters(json{
            {"beta", beta_},
            {"mu", mu_},
            {"eps", eps_}
    })
    {}

    const static std::string name() {
        return "ComplexXYModel";
    }

    typedef ComplexXYModel<SamplerCl> Model;

private:
    friend class ComplexXYModel<SamplerCl>;

    const double beta;
    const double mu;
    const double eps;
};

template<typename SamplerCl>
class ComplexXYModel : public LatticeModel< ComplexXYModel<SamplerCl> >
{
public:
    explicit ComplexXYModel(const ComplexXYModelParameters<SamplerCl> &mp_) :
        mp(mp_), sampler(SamplerCl(mp.eps)) {
    }

    // These might be the wrong distributions
    template<typename T>
    T random_state()
    {
        return normalize(sampler.template random_state<T>());
    }

    template<typename T>
    T propose_state(T site)
    {
        return normalize(sampler.template propose_state<T>(site));
    }

    /* std::complex<double> random_state()
    {
        return {random_phi_re(gen), random_phi_im(gen)};
    }

    std::complex<double> propose_state(const std::complex<double> site)
    {
        std::complex<double> state = site + std::complex<double>(propose_normal(gen), propose_normal(gen));
        return normalize(state);
    } */

    std::complex<double> normalize(std::complex<double> state)
    {
        state.real(normalize(state.real()));
        return state;
    }

    double normalize(double state)
    {
        state = std::fmod(state, 2 * M_PI);
        if (state < 0) {
            state += 2 * M_PI;
        }
        return state;
    }

    std::complex<double> get_potential(const std::complex<double> site, const std::vector<std::complex<double>*> neighbours)
    {
        double S_re = 0;
        double S_im = 0;
        for(auto i = 0; i < neighbours.size(); i += 2) {
            S_re += std::cos(site.real() - neighbours[i]->real()) * std::cosh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                    std::cos(neighbours[i+1]->real() - site.real()) * std::cosh(neighbours[i+1]->imag() - site.imag() - mp.mu * int(i == 0));
            S_im += std::sin(site.real() - neighbours[i]->real()) * std::sinh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                    std::sin(neighbours[i+1]->real() - site.real()) * std::sinh(neighbours[i+1]->imag() - site.imag() - mp.mu * int(i == 0));
        }
        return {-1.0 * mp.beta * S_re, mp.beta * S_im};
    }

    std::complex<double> get_drift_term(const std::complex<double> site, const std::vector<std::complex<double>*> neighbours)
    {
        double S_re = 0;
        double S_im = 0;
        for(auto i = 0; i < neighbours.size(); i += 2) {
            S_re += std::sin(site.real() - neighbours[i]->real()) * std::cosh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                    std::sin(site.real() - neighbours[i+1]->real()) * std::cosh(site.imag() - neighbours[i+1]->imag() + mp.mu * int(i==0));
            S_im += std::cos(site.real() - neighbours[i]->real()) * std::sinh(site.imag() - neighbours[i]->imag() - mp.mu * int(i == 0)) +
                    std::cos(site.real() - neighbours[i+1]->real()) * std::sinh(site.imag() - neighbours[i+1]->imag() + mp.mu * int(i==0));
        }
        return  mp.beta * std::complex<double>(S_re, S_im);
    }

    // [
    // For Cobrid Monte Carlo algorithms

    double get_potential(const double site, const std::vector<double*> neighbours)
    {
        std::complex<double> S = 0;
        for(auto i = 0; i < neighbours.size(); i += 2) {
            S += std::cos(std::complex<double>(site - *neighbours[i], -1.0 * mp.mu * int(i == 0))) + std::cos(std::complex<double>(*neighbours[i + 1] - site, -1.0 * mp.mu * int(i == 0)));
        }
        return (-1.0 * mp.beta * S).real(); // 0.5
    }

    double get_imag_potential(const double site, const std::vector<double*> neighbours)
    {
        std::complex<double> S = 0;
        for(auto i = 0; i < neighbours.size(); i += 2) {
            S += std::cos(std::complex<double>(site - *neighbours[i], -1.0 * mp.mu * int(i == 0))) + std::cos(std::complex<double>(*neighbours[i + 1] - site, -1.0 * mp.mu * int(i == 0)));
        }
        return (-1.0 * mp.beta * S).imag(); // 0.5
    }

    double get_drift_term(const double site, const std::vector<double*> neighbours) {
        std::complex<double> S = 0;
        for (auto i = 0; i < neighbours.size(); i += 2) {
            S += std::sin(std::complex<double>(site - *neighbours[i], mp.mu * int(i == 0))) +
                 std::sin(std::complex<double>(site - *neighbours[i + 1], -1.0 * mp.mu * int(i == 0)));
        }
        return (mp.beta * S).real();
    }

    double get_imag_drift_term(const double site, const std::vector<double*> neighbours) {
        std::complex<double> S = 0;
        for (auto i = 0; i < neighbours.size(); i += 2) {
            S += std::sin(std::complex<double>(site - *neighbours[i], mp.mu * int(i == 0))) +
                 std::sin(std::complex<double>(site - *neighbours[i + 1], -1.0 * mp.mu * int(i == 0)));
        }
        return (mp.beta * S).imag();
    }

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

private:
    const ComplexXYModelParameters<SamplerCl> &mp;

    SamplerCl sampler;
};


#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP
