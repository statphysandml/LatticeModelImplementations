//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_LANGEVIN_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_LANGEVIN_UPDATE_HPP

#include "../../mcmc_update_base.hpp"


template<typename ModelParameters, typename SamplerCl>
class NVecLangevinUpdate;


template<typename ModelParameters, typename SamplerCl>
class NVecLangevinUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit NVecLangevinUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                epsilon(get_value_by_key<double>("epsilon")), sqrt2epsilon(sqrt(2 * get_value_by_key<double>("epsilon")))
    {}

    explicit NVecLangevinUpdateParameters(
            const double epsilon_
    ) : NVecLangevinUpdateParameters(json{
            {"epsilon", epsilon_},
            {"eps", epsilon_}
    })
    {}

    static std::string name() {
        return "NVecLangevinUpdate";
    }

    typedef NVecLangevinUpdate<ModelParameters, SamplerCl> MCMCUpdate;

private:
    friend class NVecLangevinUpdate<ModelParameters, SamplerCl>;

    const double epsilon;
    const double sqrt2epsilon;
};


template<typename ModelParameters, typename SamplerCl>
class NVecLangevinUpdate : public MCMCUpdateBase< NVecLangevinUpdate<ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit NVecLangevinUpdate(const NVecLangevinUpdateParameters<ModelParameters, SamplerCl> &up_, typename ModelParameters::Model & model_)
        : MCMCUpdateBase<NVecLangevinUpdate<ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_), model(model_)
    {
        normal = std::normal_distribution<double> (0,1);
    }

    template<typename T>
    T operator() (const T site)
    {
        return update(site, model.get_drift_term(site), up.epsilon, up.sqrt2epsilon);
    }
    
private:
    const NVecLangevinUpdateParameters<ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;

    std::normal_distribution<double> normal;

    template<typename T>
    T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon)
    {
        T new_site(site);
        new_site(0) = site.orig() - epsilon * drift_term(0) + sqrt2epsilon * normal(gen);
        for(auto i = 1; i < new_site.dim(); i++)
            new_site(i) = site(i) - epsilon * drift_term(i);

        new_site.rescale(100.0);
        // new_site.remap();
        /* if(fabs(new_site.orig()) > 10)
        {
            new_site(0) = new_site.reduce();
            for(auto i = 1; i < new_site.dim(); i++)
                new_site(i) = 0.0;
        } */

        return new_site;
    }
};


#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_LANGEVIN_UPDATE_HPP
