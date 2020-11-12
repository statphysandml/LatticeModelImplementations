//
// Created by lukas on 11.11.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_VARIABLE_LANGEVIN_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_VARIABLE_LANGEVIN_UPDATE_HPP

#include "../../mcmc_update_base.hpp"


template<typename ModelParameters, typename SamplerCl>
class NVecVariableLangevinUpdate;


template<typename ModelParameters, typename SamplerCl>
class NVecVariableLangevinUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit NVecVariableLangevinUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                epsilon(get_value_by_key<double>("epsilon")), sqrt2epsilon(sqrt(2 * get_value_by_key<double>("epsilon")))
    {}

    explicit NVecVariableLangevinUpdateParameters(
            const double epsilon_
    ) : NVecVariableLangevinUpdateParameters(json{
            {"epsilon", epsilon_},
            {"eps", epsilon_}
    })
    {}

    static std::string name() {
        return "NVecVariableLangevinUpdate";
    }

    typedef NVecVariableLangevinUpdate<ModelParameters, SamplerCl> MCMCUpdate;

private:
    friend class NVecVariableLangevinUpdate<ModelParameters, SamplerCl>;

    const double epsilon;
    const double sqrt2epsilon;
};


template<typename ModelParameters, typename SamplerCl>
class NVecVariableLangevinUpdate : public MCMCUpdateBase< NVecVariableLangevinUpdate<ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit NVecVariableLangevinUpdate(const NVecVariableLangevinUpdateParameters<ModelParameters, SamplerCl> &up_, typename ModelParameters::Model & model_)
            : MCMCUpdateBase<NVecVariableLangevinUpdate<ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_), model(model_)
    {
        normal = std::normal_distribution<double> (0,1);
    }

    template<typename T>
    T operator() (const T site)
    {
        return update(site, model.get_drift_term(site), up.epsilon, up.sqrt2epsilon);
    }

private:
    const NVecVariableLangevinUpdateParameters<ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;

    std::normal_distribution<double> normal;

    template<typename T>
    T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon)
    {
        model.update_ft();

        T new_site(site);
        new_site(0) = site.orig() - epsilon * drift_term(0) + sqrt2epsilon * normal(gen);
        for(uint i = 1; i < new_site.dim(); i++)
            new_site(i) = site(i) - epsilon * drift_term(i);

        new_site.rescale(10.0);

        return new_site;
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_VARIABLE_LANGEVIN_UPDATE_HPP
