//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMIAL_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"

#include "../site_model.hpp"


template<typename SamplerCl, typename NVecSplitting>
class NVecPolynomialModel;

template<typename SamplerCl, typename NVecSplitting>
class NVecPolynomialModelParameters : public SiteModelParameters {
public:
    explicit NVecPolynomialModelParameters<SamplerCl, NVecSplitting>(const json params_) : SiteModelParameters(params_),
                                                                            sigma(get_value_by_key< double >("sigma")),
                                                                            lambda(get_value_by_key< double >("lambda")),
                                                                            h(get_value_by_key< double >("h"))
    {}

    explicit NVecPolynomialModelParameters<SamplerCl, NVecSplitting>(double beta_, double mu_, double eps_) : NVecPolynomialModelParameters(json {})
    {}

    const static std::string name() {
        return "NVecPolynomialModel";
    }

    typedef NVecPolynomialModel<SamplerCl, NVecSplitting> Model;

private:
    friend class NVecPolynomialModel<SamplerCl, NVecSplitting>;

    const double sigma;
    const double lambda;
    const double h;
};


template<typename SamplerCl, typename NVecSplitting>
class NVecPolynomialModel : public SiteModel< NVecPolynomialModel<SamplerCl, NVecSplitting> >
{
public:
    explicit NVecPolynomialModel(const NVecPolynomialModelParameters<SamplerCl, NVecSplitting> &mp_) :
        mp(mp_), nvec_splitting(NVecSplitting(mp.lambda, mp.sigma, mp.h))
    {}

    using SiteModel< NVecPolynomialModel<SamplerCl, NVecSplitting> >::random_state;

    /* typename NVecSplitting::Ttype random_state()
    {
        return typename NVecSplitting::Ttype(sampler.template random_state<typename NVecSplitting::Ttype>());
    }

    typename NVecSplitting::Ttype propose_state(typename NVecSplitting::Ttype site)
    {
        return typename NVecSplitting::Ttype(sampler.template propose_state<typename NVecSplitting::Ttype>(), site);
    } */

    typename NVecSplitting::Ttype get_potential(const typename NVecSplitting::Ttype site) const
    {
        return typename NVecSplitting::Ttype(nvec_splitting.get_potential(site));
    };

    typename NVecSplitting::Ttype get_drift_term(const typename NVecSplitting::Ttype site) const
    {
        return typename NVecSplitting::Ttype(nvec_splitting.get_drift_term(site));
    };

private:
    const NVecPolynomialModelParameters<SamplerCl, NVecSplitting> &mp;
    const NVecSplitting nvec_splitting;
};

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMIAL_MODEL_HPP
