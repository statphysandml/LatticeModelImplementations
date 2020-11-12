//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMIAL_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"

#include "../site_model.hpp"


template<typename NVecSplitting>
class NVecPolynomialModel;

template<typename NVecSplitting>
class NVecPolynomialModelParameters : public SiteModelParameters {
public:
    explicit NVecPolynomialModelParameters<NVecSplitting>(const json params_) : SiteModelParameters(params_),
                                                                            sigma(get_value_by_key< double >("sigma")),
                                                                            lambda(get_value_by_key< double >("lambda")),
                                                                            h(get_value_by_key< double >("h")),
                                                                            drift_term_in_python_code(NVecSplitting::get_drift_term_in_python_code())
    {}

    explicit NVecPolynomialModelParameters<NVecSplitting>(
        double lambda_, double sigma_, double h_) : NVecPolynomialModelParameters(json {
            {"lambda", lambda_},
            {"sigma", sigma_},
            {"h", h_},
            {"drift_term_in_python_code", NVecSplitting::get_drift_term_in_python_code()}
    })
    {}

    const static std::string name() {
        return "NVecPolynomialModel" + NVecSplitting::name();
    }

    typedef NVecPolynomialModel<NVecSplitting> Model;

private:
    friend class NVecPolynomialModel<NVecSplitting>;

    const double sigma;
    const double lambda;
    const double h;
    const std::string drift_term_in_python_code;
};


template<typename NVecSplitting>
class NVecPolynomialModel : public SiteModel< NVecPolynomialModel<NVecSplitting> >
{
public:
    explicit NVecPolynomialModel(const NVecPolynomialModelParameters<NVecSplitting> &mp_) :
        mp(mp_), nvec_splitting(NVecSplitting(mp.lambda, mp.sigma, mp.h))
    {}

    template<typename T>
    T get_potential(const T site) const
    {
        return T(nvec_splitting.get_potential(site));
    };

    template<typename T>
    T get_drift_term(const T site) const
    {
        return T(nvec_splitting.get_drift_term(site));
    };

private:
    const NVecPolynomialModelParameters<NVecSplitting> &mp;
    const NVecSplitting nvec_splitting;
};

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_POLYNOMIAL_MODEL_HPP
