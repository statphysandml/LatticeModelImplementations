//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_SU_TWO_MODEL_HPP
#define MAIN_SU_TWO_MODEL_HPP


#include "../link_lattice_model.hpp"

// #include "../su2.hpp"

class SUTwoModel;

class SUTwoModelParameters : public LinkLatticeModelParameters {
public:
    explicit SUTwoModelParameters(const json params_) : LinkLatticeModelParameters(params_),
                                                       beta(get_value_by_key<double>("beta")),
                                                       eps(get_value_by_key<double>("eps"))
    {}

    explicit SUTwoModelParameters(double beta_, double eps_) : SUTwoModelParameters(json{
            {"beta", beta_},
            {"eps", eps_}
    })
    {}

    const static std::string name() {
        return "SUTwoModel";
    }

    static uint N() {
        return 2;
    };

    typedef SUTwoModel Model;

private:
    friend class SUTwoModel;

    const double beta;
    const double eps;
};


class SUTwoModel : public LinkLatticeModel< SUTwoModel >
{
public:
    explicit SUTwoModel(const SUTwoModelParameters &mp_) : mp(mp_) {}

    template< typename T>
    T random_state()
    {
        return T("random");
    }

    template<typename T>
    T propose_state(T site)
    {
        return site * T(mp.eps);
    }

    template<typename T, typename T2=double_t>
    T2 get_potential(const T site, const std::vector<T*> neighbours)
    {
        T A = calc_A(neighbours);
        // (lat.dim()-1)*2-std::real(1.0/N*(lat(n*lat.dim()+mu)*calc_A(lat, n, mu, both_orientations)).trace())
        return mp.beta/SUTwoModelParameters::N() * (neighbours.size() / 3 * SUTwoModelParameters::N() - std::real((site*A).trace())); // Additional 1/N missing in second term??
    }

private:
    const SUTwoModelParameters &mp;
};


#endif //MAIN_SU_TWO_MODEL_HPP
