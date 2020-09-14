//
// Created by lukas on 05.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_IMAGINARY_GAUSSIAN_DISTRIBUTION_HPP
#define LATTICEMODELIMPLEMENTATIONS_IMAGINARY_GAUSSIAN_DISTRIBUTION_HPP


#include "../site/site.hpp"
#include "../update_dynamics/site_update_dynamics/site_update_formalisms/simple_update.hpp"
#include "../mcmc_update/local_mcmc_update/local_update_formalisms/complex_langevin_update.hpp"
#include "../site/site_models/complex_scalar_gaussian_model.hpp"


class imaginary_gaussian_distribution
{
public:
    imaginary_gaussian_distribution(std::complex<double> a_, std::complex<double> b_, uint n_autocorrelation_, uint n_initialization, const double epsilon) :
        a(a_), b(b_), n_autocorrelation(n_autocorrelation_),
        complex_gaussian_site_system_parameters(SiteParameters< std::complex<double>, ModelParams, ComplexLangevinUpdateParameters<ModelParams>, SiteSimpleUpdateParameters>(
                json {
                    {"measures", {"Mean", "ComplexConfig"}},
                    {ModelParams::param_file_name(), {{"a", a_}, {"b", b_}}},
                    {ComplexLangevinUpdateParameters<ModelParams>::param_file_name(), {{"epsilon", epsilon}}},
                    {SiteSimpleUpdateParameters::param_file_name(), {}}},
                    "None")),
        complex_gaussian_site_system(complex_gaussian_site_system_parameters)
    {
        complex_gaussian_site_system.update_step(n_initialization);
    }

    /* template< class Generator >
    std::complex<double> operator() (Generator& g )
    {
        complex_gaussian_site_system.update_step(n_autocorrelation);
        return complex_gaussian_site_system.get_system_representation();
    } */

    template< class Generator >
    double operator()( Generator& g )
    {
        complex_gaussian_site_system.update_step(n_autocorrelation);
        return complex_gaussian_site_system.get_system_representation().real();
    }

    const static std::string name()
    {
        return "imaginary_gaussian_distribution";
    }

private:
    typedef ComplexScalarGaussianModelParameters<GaussianSampler> ModelParams;
    const SiteParameters< std::complex<double>, ModelParams, ComplexLangevinUpdateParameters<ModelParams>, SiteSimpleUpdateParameters> complex_gaussian_site_system_parameters;
    SiteSystem< std::complex<double>, ModelParams, ComplexLangevinUpdateParameters<ModelParams>, SiteSimpleUpdateParameters> complex_gaussian_site_system;

    const std::complex<double> a;
    const std::complex<double> b;
    const uint n_autocorrelation;
};

#endif //LATTICEMODELIMPLEMENTATIONS_IMAGINARY_GAUSSIAN_DISTRIBUTION_HPP
