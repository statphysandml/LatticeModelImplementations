//
// Created by lukas on 10.10.19.
//

#ifndef MAIN_LANGEVIN_UPDATE_HPP
#define MAIN_LANGEVIN_UPDATE_HPP

#include "../local_update.hpp"


template<typename ModelParameters>
class LangevinUpdate;


template<typename ModelParameters>
class LangevinUpdateParameters : public LocalUpdateFormalismParameters {
public:
    explicit LangevinUpdateParameters(const json params_) : LocalUpdateFormalismParameters(params_),
        epsilon(get_value_by_key<double>("epsilon")),
        sqrt2epsilon(sqrt(2 * get_value_by_key<double>("epsilon")))
    {}

    explicit LangevinUpdateParameters(
            const double epsilon_
    ) : LangevinUpdateParameters(json{
            {"epsilon", epsilon_}
    })
    {}

    static std::string name() {
        return "LangevinUpdate";
    }

    typedef LangevinUpdate<ModelParameters> UpdateFormalism;

protected:
    friend class LangevinUpdate<ModelParameters>;

    const double epsilon;
    const double sqrt2epsilon;
};


template<typename ModelParameters>
class LangevinUpdate : public LocalUpdateFormalism< LangevinUpdate<ModelParameters> >
{
public:
    explicit LangevinUpdate(const LangevinUpdateParameters<ModelParameters> &up_, typename ModelParameters::Model & model_) : up(up_), model(model_)
    {
        normal = std::normal_distribution<double> (0,1);
    }

    template<typename Lattice>
    void initialize_local_update(const Lattice &lattice)
    {}

    template<typename T>
    T operator() (const T site)
    {
        const T drift_term = model.get_drift_term(site);
        return model.normalize(site - up.epsilon * drift_term + up.sqrt2epsilon * normal(gen));
    }

    template<typename T>
    T operator() (const T site, const double KMax, const double KExpectation)
    {
        double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);

        const T drift_term = model.get_drift_term(site);
        return model.normalize(site - epsilon * drift_term + std::sqrt(2 * epsilon) * normal(gen));
    }

    template<typename T>
    T operator() (const T site, const std::vector< T* > neighbours)
    {
        const T drift_term = model.get_drift_term(site, neighbours);
        return model.normalize(site - up.epsilon * drift_term + up.sqrt2epsilon * normal(gen));
    }

    template<typename T>
    T operator() (const T site, const std::vector< T* > neighbours, const double KMax, const double KExpectation)
    {
        double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);

        const T drift_term = model.get_drift_term(site, neighbours);
        return model.normalize(site - epsilon * drift_term + std::sqrt(2 * epsilon) * normal(gen));
    }

protected:
    const LangevinUpdateParameters<ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::normal_distribution<double> normal;
};

#endif //MAIN_LANGEVIN_UPDATE_HPP
