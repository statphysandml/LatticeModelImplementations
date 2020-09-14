//
// Created by lukas on 31.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP


#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"


namespace xymodel
{
    template<typename SB>
    struct MeasureXYMagnetization: public common_measures::MeasurePolicy< SB > {
    public:
        std::string measure(const SB &system) override {
            auto sum_cos = decltype(system[0]){0};
            auto sum_sin = decltype(system[0]){0};

            for(auto i = 0; i < system.size(); i++) {
                sum_cos += std::cos(system[i]);
                sum_sin += std::sin(system[i]);
            }

            return std::to_string(std::pow(sum_cos / system.size(), 2) + std::pow(sum_sin / system.size(), 2));
        }

        std::string name()
        {
            return "XYMagnetization";
        }
    };
}


template<typename SamplerCl>
class XYModel;

template<typename SamplerCl>
class XYModelParameters : public LatticeModelParameters {
public:
    explicit XYModelParameters(const json params_) : LatticeModelParameters(params_),
                                                     beta(get_value_by_key<double>("beta", 0.5)),
                                                     J(get_value_by_key<double>("J", 1.0)),
                                                     h(get_value_by_key<double>("h", 0.0)),
                                                     eps(get_value_by_key<double>("eps", 0.1))
    {}

    explicit XYModelParameters(double beta_, double J_, double h_, double eps_) : XYModelParameters(json{
            {"beta", beta_},
            {"J", J_},
            {"h", h_},
            {"eps", eps_}
    })
    {}

    const static std::string name() {
        return "XYModel";
    }

    typedef XYModel<SamplerCl> Model;

private:
    friend class XYModel<SamplerCl>;

    const double beta;
    const double J;
    const double h;
    const double eps;
};


template<typename SamplerCl>
class XYModel : public LatticeModel< XYModel<SamplerCl> >
{
public:
    explicit XYModel(const XYModelParameters<SamplerCl> &mp_) : mp(mp_), sampler(SamplerCl(mp.eps)) {}

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

    template<typename T>
    T normalize(T state)
    {
        state = std::fmod(state, 2 * M_PI);
        if (state < 0) {
            state += 2 * M_PI;
        }
        return state;
    }

    template<typename T>
    T get_potential(const T site, const std::vector<T*> neighbours)
    {
        double S = 0;
        for(auto i = 0; i < neighbours.size(); i += 2) {
            S += mp.J * std::cos(site - *neighbours[i]) + mp.J * std::cos(*neighbours[i + 1] - site) + mp.h * std::cos(site);
        }
        return -1.0 * mp.beta * S; // 0.5
    }

    template<typename T>
    T get_drift_term(const T site, const std::vector<T*> neighbours)
    {
        double S = 0;
        for(auto i = 0; i < neighbours.size(); i += 2) {
            S += mp.J * std::sin(site - *neighbours[i]) +
                    mp.J * std::sin(site - *neighbours[i+1]) + mp.h * std::cos(site);
        }
        return  mp.beta * S;
    }

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

    template<typename SB>
    std::vector< std::unique_ptr<common_measures::MeasurePolicy<SB>> > generate_model_measures(const json& measure_names)
    {
        std::vector< std::unique_ptr<common_measures::MeasurePolicy<SB>> > measures {};
        for (auto& measure_name :  measure_names)
            if(measure_name == "XYMagnetization")
                measures.push_back(std::make_unique<xymodel::MeasureXYMagnetization<SB>>());
        return measures;
    }


private:
    const XYModelParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

#endif //LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP
