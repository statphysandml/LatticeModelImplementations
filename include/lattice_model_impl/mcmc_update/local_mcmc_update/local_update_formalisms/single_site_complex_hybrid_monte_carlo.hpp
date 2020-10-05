//
// Created by lukas on 24.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_HYBRID_MONTE_CARLO_HPP

#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"
#include "boost/proto/functional/std/utility.hpp"

#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"

#include "../../mcmc_update_base.hpp"

#ifdef THRUST
#include "../../../thrust/thrust_complex_gaussian_distribution.hpp"
#endif

#include "../../../util/distribution/imaginary_gaussian_distribution.hpp"
#include "../../../util/distribution/complex_gaussian_distribution_from_file.hpp"


template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexHybridMonteCarloUpdate;

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexHybridMonteCarloUpdateParameters;

template<typename ModelParameters, typename MomentumDistribution>
struct momentum_distribution
{
    template<template<typename, typename> class UpdateFormalism>
    explicit momentum_distribution(const UpdateFormalism<ModelParameters, MomentumDistribution> &up_) :
            sampler(std::complex<double>{up_.dist_sigma_real, up_.dist_sigma_imag}, up_.dist_sampler_n_autocorrelation, up_.dist_sampler_n_initialization, up_.dist_sampler_epsilon, 10100000)
    {}

    template< class Generator >
    double operator()( Generator& g )
    {
        return sampler(g);
    }

    template< class Generator >
    std::complex<double> sample_complex(Generator& g)
    {
        return sampler.sample_complex();
    }

    MomentumDistribution sampler;
};

template<typename ModelParameters, typename SamplerCl>
struct momentum_distribution<ModelParameters, std::normal_distribution<double>>
{
    template<template<typename, typename> class UpdateFormalism>
    explicit momentum_distribution(const UpdateFormalism<ModelParameters, std::normal_distribution<double>> &up_) :
            sampler(0.0, 1.0)
    {}

    template< class Generator >
    double operator()( Generator& g )
    {
        return sampler(g);
    }

    template< class Generator >
    std::complex<double> sample_complex(Generator& g)
    {
        return {sampler(g), 0.0};
    }

    std::normal_distribution<double> sampler;
};

template<typename ModelParameters, typename SamplerCl>
struct momentum_distribution<ModelParameters, complex_gaussian_distribution_from_file>
{
    template<template<typename, typename> class UpdateFormalism>
    explicit momentum_distribution(const UpdateFormalism<ModelParameters, complex_gaussian_distribution_from_file> &up_) :
            sampler(up_.get_sigma(), 21000000, "ThrustGaussianDistribution")
    {}

    template< class Generator >
    double operator()( Generator& g )
    {
        return sampler(g);
    }

    template< class Generator >
    std::complex<double> sample_complex(Generator& g)
    {
        return sampler.sample_complex();
    }

    complex_gaussian_distribution_from_file sampler;
};

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexHybridMonteCarloUpdateParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteComplexHybridMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                              dt(get_value_by_key<double>("dt", 0.01)),
                                                                              n(get_value_by_key<int>("n", 20)),
                                                                              dist_sigma_real(get_value_by_key<double>("dist_sigma_real", 0.001)),
                                                                              dist_sigma_imag(get_value_by_key<double>("dist_sigma_imag", 1.0)),
                                                                              dist_sampler_n_autocorrelation(get_value_by_key<uint>("dist_sampler_n_autocorrelation", 10000)),
                                                                              dist_sampler_n_initialization(get_value_by_key<uint>("dist_sampler_n_initialization", 10000)),
                                                                              dist_sampler_epsilon(get_value_by_key<double>("dist_sampler_epsilon", 0.0001)),
                                                                              reweight_momentum(get_value_by_key<std::string>("reweight_momentum", "true"))


    {}

    explicit SingleSiteComplexHybridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : SingleSiteComplexHybridMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteComplexHybridMonteCarloUpdate";
    }

    std::complex<double> get_sigma() const
    {
        return {dist_sigma_real, dist_sigma_imag};
    }

    typedef SingleSiteComplexHybridMonteCarloUpdate<ModelParameters, MomentumDistribution> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteComplexHybridMonteCarloUpdate<ModelParameters, MomentumDistribution>;
    friend struct momentum_distribution<ModelParameters, MomentumDistribution>;

    const double dt;
    const int n;
    const double dist_sigma_real;
    const double dist_sigma_imag;
    const uint dist_sampler_n_autocorrelation;
    const uint dist_sampler_n_initialization;
    const double dist_sampler_epsilon;
    const std::string reweight_momentum;
};

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexHybridMonteCarloUpdate : public MCMCUpdateBase< SingleSiteComplexHybridMonteCarloUpdate<ModelParameters, MomentumDistribution> >
{
public:
    explicit SingleSiteComplexHybridMonteCarloUpdate(const SingleSiteComplexHybridMonteCarloUpdateParameters<ModelParameters, MomentumDistribution> &up_, typename ModelParameters::Model & model_) :
            up(up_), model(model_), gaussian_dist(up)
    {
        proposal_site = std::vector<double>{0.0, 0.0, 0.0, 0.0};
        current_site = std::vector<double>{0.0, 0.0, 0.0, 0.0};
    }

    std::complex<double> operator() (const std::complex<double> site)
    {
        std::cout << "Site " << site.real() << " + i" << site.imag() << std::endl;

        // X
        proposal_site[0] = site.real();
        proposal_site[1] = site.imag();

        // Sample momenta
        auto sample = gaussian_dist.sample_complex(gen);
        proposal_site[2] = sample.real();
        proposal_site[3] = sample.imag();

        current_site = proposal_site;

        hamilton_system_H_imag_partially_constant hamiltonian_system(model, up.dist_sigma_real, up.dist_sigma_imag);

        auto energy = hamiltonian_system.get_energy(proposal_site);

        boost::numeric::odeint::integrate_n_steps(stepper_type(),
                                                  hamiltonian_system,
                                                  proposal_site,
                                                  0.0, up.dt, up.n);

        auto proposal_energy = hamiltonian_system.get_energy(proposal_site);

        // std::cout << proposal_energy.real() - energy.real() << " + i" << proposal_energy.imag() - energy.imag() << std::endl;

        return hamiltonian_system.accept_reject_step(proposal_site, current_site, proposal_energy, energy);
    }

    typedef boost::numeric::odeint::runge_kutta4< std::vector<double> > stepper_type;
private:

    struct hamilton_system_H_imag_constant
    {
        hamilton_system_H_imag_constant(typename ModelParameters::Model & model_, const double real_sigma_, const double imag_sigma_) : hmc_model(model_), real_sigma(real_sigma_), imag_sigma(imag_sigma_)
        {
            rand = std::uniform_real_distribution<double> (0,1);
        }

        void operator()( const std::vector<double> &x , std::vector<double> &dxdt, const double /* t */  ) const {
            auto drift = hmc_model.get_drift_term(std::complex<double> {x[0], x[1]});
            dxdt[2] = -1.0 * drift.imag();
            dxdt[3] = -1.0 * drift.real();
            dxdt[0] = real_sigma * x[3] + imag_sigma * x[2];
            dxdt[1] = real_sigma * x[2] - imag_sigma * x[3];
        }

        std::complex<double> get_energy(std::vector<double> &x) const
        {
            return hmc_model.get_potential(std::complex<double>{x[0], x[1]}) + (0.5 * std::complex<double> (real_sigma, imag_sigma) * std::pow(std::complex<double> {x[2], x[3]}, 2.0));
        }

        std::complex<double> accept_reject_step(std::vector<double> new_site, std::vector<double> site, std::complex<double> proposal_energy, std::complex<double> energy)
        {
            if(rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy.real() - energy.real())))) //  and std::abs(proposal_energy.imag() - current_energy.imag()) < up.max_delta_s_imag)
                return {new_site[0], new_site[1]};
            else
                return {site[0], site[1]};
        }

        typename ModelParameters::Model & hmc_model;
        const double real_sigma;
        const double imag_sigma;
        std::uniform_real_distribution<double> rand;
    };

    struct hamilton_system_H_constant
    {
        hamilton_system_H_constant(typename ModelParameters::Model & model_, const double real_sigma_, const double imag_sigma_) : hmc_model(model_), real_sigma(real_sigma_), imag_sigma(imag_sigma_)
        {}

        void operator()( const std::vector<double> &x , std::vector<double> &dxdt, const double /* t */  ) const {
            auto drift = hmc_model.get_drift_term(std::complex<double> {x[0], x[1]});
            // ist genau andersherum wie H_imag constant
            dxdt[2] = -1.0 * drift.real();
            dxdt[3] = -1.0 * drift.imag();
            dxdt[0] = real_sigma * x[2] - imag_sigma * x[3];
            dxdt[1] = real_sigma * x[3] + imag_sigma * x[2];
        }

        std::complex<double> get_energy(std::vector<double> &x) const
        {
            return hmc_model.get_potential(std::complex<double>{x[0], x[1]}) + (0.5 * std::complex<double> (real_sigma, imag_sigma) * std::pow(std::complex<double> {x[2], x[3]}, 2.0));
        }

        std::complex<double> accept_reject_step(std::vector<double> new_site, std::vector<double> site, std::complex<double> proposal_energy, std::complex<double> energy)
        {
            return {new_site[0], new_site[1]};
        }

        typename ModelParameters::Model & hmc_model;
        const double real_sigma;
        const double imag_sigma;
    };

    struct hamilton_system_H_imag_partially_constant
    {
        hamilton_system_H_imag_partially_constant(typename ModelParameters::Model & model_, const double real_sigma_, const double imag_sigma_) : hmc_model(model_), real_sigma(real_sigma_), imag_sigma(imag_sigma_)
        {
            rand = std::uniform_real_distribution<double> (0,1);
        }

        void operator()( const std::vector<double> &x , std::vector<double> &dxdt, const double /* t */  ) const {
            auto drift = hmc_model.get_complex_hybrid_drift_term(std::complex<double> {x[0], x[1]});
            dxdt[2] = -1.0 * drift.real();
            dxdt[3] = -1.0 * drift.imag();
            dxdt[0] = imag_sigma * x[2];
            dxdt[1] = imag_sigma * x[3];
        }

        std::complex<double> get_energy(std::vector<double> &x) const
        {
            return hmc_model.get_complex_hybrid_potential(std::complex<double>{x[0], x[1]}) + (0.5 * std::complex<double> (imag_sigma, 0.0) * std::pow(std::complex<double> {x[2], x[3]}, 2.0));
        }

        std::complex<double> accept_reject_step(std::vector<double> new_site, std::vector<double> site, std::complex<double> proposal_energy, std::complex<double> energy)
        {
            auto current_energy = hmc_model.get_potential({site[0], site[1]});
            auto new_energy = hmc_model.get_potential({new_site[0], site[1]});

            if (rand(gen) < std::min(
                    1.0, std::exp(-1.0 * (new_energy - current_energy).real() - real_sigma / 2.0 * (std::pow(std::complex<double>{new_site[2], new_site[3]}, 2.0) - std::pow(std::complex<double>{site[2], site[3]}, 2.0)).real())))
            {
                new_site[1] = site[1] - 0.01 * (new_energy - current_energy).imag() / (new_site[0] - site[0]);
                return {new_site[0], new_site[1]};
            }
            else
                return {site[0], site[1]};
        }

        typename ModelParameters::Model & hmc_model;
        const double real_sigma;
        const double imag_sigma;
        std::uniform_real_distribution<double> rand;
    };

    const SingleSiteComplexHybridMonteCarloUpdateParameters<ModelParameters, MomentumDistribution> & up;
    typename ModelParameters::Model & model;

    std::vector<double> proposal_site;
    std::vector<double> current_site;

    momentum_distribution<ModelParameters, MomentumDistribution> gaussian_dist;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_HYBRID_MONTE_CARLO_HPP
