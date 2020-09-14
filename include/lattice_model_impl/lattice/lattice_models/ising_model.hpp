//
// Created by lukas on 17.10.19.
//

#ifndef MAIN_ISING_MODEL_HPP
#define MAIN_ISING_MODEL_HPP

#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"

class IsingModel;

class IsingModelParameters : public LatticeModelParameters {
public:
    explicit IsingModelParameters(const json params_) : LatticeModelParameters(params_),
                                                        beta(get_value_by_key<double>("beta", 0.4)),
                                                        J(get_value_by_key<std::complex<double>>("J", 1.0)),
                                                        h(get_value_by_key<std::complex<double>>("h", 0.0))
    {

    }

    explicit IsingModelParameters(double beta_, std::complex<double> J_, std::complex<double> h_) : IsingModelParameters(json{
            {"beta", beta_},
            {"J", J_},
            {"h", h_}
    })
    {}

    const static std::string name() {
        return "IsingModel";
    }

    typedef IsingModel Model;

private:
    friend class IsingModel;

    const double beta;
    const std::complex<double> J;
    const std::complex<double> h;
};


class IsingModel : public LatticeModel< IsingModel >
{
public:
    explicit IsingModel(const IsingModelParameters &mp_) : mp(mp_) {
        uniint = std::uniform_int_distribution<int>(0, 1);
    }

    template< typename T>
    T random_state()
    {
        return 2 * uniint(gen) - 1;
    }

    template<typename T>
    T propose_state(T site)
    {
        return -1 * site;
    }

    template<typename T, typename T2=double_t>
    T2 get_potential(const T site, const std::vector<T*> neighbours)
    {
        double coupling = 0;
        for(auto i = 0; i < neighbours.size(); i++) {
            coupling += *neighbours[i];
        }
        return  -1.0 * mp.beta * site * (mp.J.real() * coupling + mp.h.real()); // 0.5
    }

private:
    const IsingModelParameters &mp;
    std::uniform_int_distribution<int> uniint;
};


#endif //MAIN_ISING_MODEL_HPP
