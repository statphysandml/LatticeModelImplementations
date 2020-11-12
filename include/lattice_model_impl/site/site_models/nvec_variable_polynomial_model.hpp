//
// Created by lukas on 11.11.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_VARIABLE_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_VARIABLE_POLYNOMIAL_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"

#include "../site_model.hpp"


template<typename NVecSplitting>
class NVecVariablePolynomialModel;

template<typename NVecSplitting>
class NVecVariablePolynomialModelParameters : public SiteModelParameters {
public:
    explicit NVecVariablePolynomialModelParameters<NVecSplitting>(const json params_) : SiteModelParameters(params_),
                                                                                sigma(get_value_by_key< double >("sigma")),
                                                                                lambda(get_value_by_key< double >("lambda")),
                                                                                h(get_value_by_key< double >("h")),
                                                                                drift_term_in_python_code(NVecSplitting::get_drift_term_in_python_code())
    {}

    explicit NVecVariablePolynomialModelParameters<NVecSplitting>(
            double lambda_, double sigma_, double h_) : NVecVariablePolynomialModelParameters(json {
            {"lambda", lambda_},
            {"sigma", sigma_},
            {"h", h_},
            {"drift_term_in_python_code", NVecSplitting::get_drift_term_in_python_code()}
    })
    {}

    const static std::string name() {
        return "NVecVariablePolynomialModel" + NVecSplitting::name();
    }

    typedef NVecVariablePolynomialModel<NVecSplitting> Model;

private:
    friend class NVecVariablePolynomialModel<NVecSplitting>;

    const double sigma;
    const double lambda;
    const double h;
    const std::string drift_term_in_python_code;
};


template<typename NVecSplitting>
class NVecVariablePolynomialModel : public SiteModel< NVecVariablePolynomialModel<NVecSplitting> >
{
public:
    explicit NVecVariablePolynomialModel(const NVecVariablePolynomialModelParameters<NVecSplitting> &mp_) :
            mp(mp_), nvec_splitting(NVecSplitting(mp.lambda, mp.sigma, mp.h))
    {}

    void update_ft()
    {
        nvec_splitting.update_ft();
    }

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
    const NVecVariablePolynomialModelParameters<NVecSplitting> &mp;
    NVecSplitting nvec_splitting;
};

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_VARIABLE_POLYNOMIAL_MODEL_HPP
