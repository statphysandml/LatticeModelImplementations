//
// Created by lukas on 04.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP

#include "../site_model.hpp"
#include "mcmc_simulation/util/random.hpp"


template<typename SamplerCl>
class ComplexScalarGaussianModel;

template<typename SamplerCl>
class ComplexScalarGaussianModelParameters : public SiteModelParameters {
public:
    explicit ComplexScalarGaussianModelParameters(const json params_) : SiteModelParameters(params_),
                                                                        eps(get_value_by_key<double>("eps", 0.1)),
                                                                        a(get_value_by_key< std::complex<double> >("a")),
                                                                        b(get_value_by_key< std::complex<double> >("b"))
    {}

    explicit ComplexScalarGaussianModelParameters(std::complex<double> a_, std::complex<double> b_, double eps_=0.1) : ComplexScalarGaussianModelParameters(json {
            {"a", a_},
            {"b", b_},
            {"eps", eps_}
    })
    {}

    static std::string name() {
        return "ComplexScalarGaussianModel";
    }

    typedef ComplexScalarGaussianModel<SamplerCl> Model;

private:
    friend class ComplexScalarGaussianModel<SamplerCl>;

    const double eps;
    const std::complex<double> a;
    const std::complex<double> b;
};

template<typename SamplerCl>
class ComplexScalarGaussianModel : public SiteModel< ComplexScalarGaussianModel<SamplerCl> >
{
public:
    explicit ComplexScalarGaussianModel(const ComplexScalarGaussianModelParameters<SamplerCl> &mp_) : mp(mp_), sampler(SamplerCl(mp.eps)) {}

    template<typename T>
    T random_state()
    {
        return sampler.template random_state<T>();
    }

    template<typename T>
    T propose_state(T site)
    {
        return sampler.template propose_state<T>(site);
    }

    std::complex<double> get_drift_term(const std::complex<double> site)
    {
        return  {mp.a.real() * site.real() - mp.a.imag() * site.imag() - mp.b.real(),
                 mp.a.real() * site.imag() + mp.a.imag() * site.real() - mp.b.imag()};
    }

    std::complex<double> get_potential(const std::complex<double> site)
    {
        return {0.5 * mp.a.real() * std::pow(site.real(), 2) - 0.5 * mp.a.real() * std::pow(site.imag(), 2) - mp.a.imag() * site.real() * site.imag()  - mp.b.real() * site.real() + mp.b.imag() * site.imag(),
                0.5 * mp.a.imag() * std::pow(site.real(), 2) - 0.5 * mp.a.imag() * std::pow(site.imag(), 2) + mp.a.real() * site.real() * site.imag()  - mp.b.real() * site.imag() - mp.b.imag() * site.real()};
    }

    const SamplerCl& get_sampler() const
    {
        return sampler;
    }

private:
    const ComplexScalarGaussianModelParameters<SamplerCl> &mp;

    SamplerCl sampler;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_SCALAR_GAUSSIAN_MODEL_HPP
