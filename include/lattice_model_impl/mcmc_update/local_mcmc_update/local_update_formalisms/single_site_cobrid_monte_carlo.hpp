//
// Created by lukas on 07.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COBRID_MONTE_CARLO_HPP


#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"
#include "boost/proto/functional/std/utility.hpp"

#include "../../mcmc_update_base.hpp"

#ifdef THRUST
#include "../../../thrust/thrust_complex_gaussian_distribution.hpp"
#endif

#include "../../../util/distribution/imaginary_gaussian_distribution.hpp"
#include "../../../util/distribution/complex_gaussian_distribution_from_file.hpp"


// Pure real approach for the outcome of x - fails

template<typename T, typename ModelParameters, typename MomentumDistribution>
class SingleSiteCobridMonteCarloUpdate;

template<typename T, typename ModelParameters, typename MomentumDistribution>
class SingleSiteCobridMonteCarloUpdateParameters;

template<typename T, typename ModelParameters, typename MomentumDistribution>
std::unique_ptr<MomentumDistribution> initialize_momentum_distribution(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, MomentumDistribution> &up_);

template<typename T, typename ModelParameters>
std::unique_ptr<std::normal_distribution<double>> initialize_momentum_distribution(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, std::normal_distribution<double> > &up_);

template<typename T, typename ModelParameters>
std::unique_ptr<complex_gaussian_distribution_from_file> initialize_momentum_distribution(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, complex_gaussian_distribution_from_file > &up_);

template<typename T, typename ModelParameters, typename MomentumDistribution>
class SingleSiteCobridMonteCarloUpdateParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteCobridMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                    dt(get_value_by_key<double>("dt", 0.01)),
                                                                    n(get_value_by_key<int>("n", 20)),
                                                                    dist_sigma_real(get_value_by_key<double>("dist_sigma_real", 0.001)),
                                                                    dist_sigma_imag(get_value_by_key<double>("dist_sigma_imag", 1.0)),
                                                                    dist_sampler_n_autocorrelation(get_value_by_key<uint>("dist_sampler_n_autocorrelation", 10000)),
                                                                    dist_sampler_n_initialization(get_value_by_key<uint>("dist_sampler_n_initialization", 10000)),
                                                                    dist_sampler_epsilon(get_value_by_key<double>("dist_sampler_epsilon", 0.0001)),
                                                                    reweight_momentum(get_value_by_key<std::string>("reweight_momentum", "true"))


    {}

    explicit SingleSiteCobridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : SingleSiteCobridMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteCobridMonteCarloUpdate";
    }

    std::complex<double> get_sigma() const
    {
        return {dist_sigma_real, dist_sigma_imag};
    }

    typedef SingleSiteCobridMonteCarloUpdate<T, ModelParameters, MomentumDistribution> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteCobridMonteCarloUpdate<T, ModelParameters, MomentumDistribution>;
    friend std::unique_ptr<MomentumDistribution> initialize_momentum_distribution<T, ModelParameters, MomentumDistribution>(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, MomentumDistribution> &up_);

    const double dt;
    const int n;
    const double dist_sigma_real;
    const double dist_sigma_imag;
    const uint dist_sampler_n_autocorrelation;
    const uint dist_sampler_n_initialization;
    const double dist_sampler_epsilon;
    const std::string reweight_momentum;
};

template<typename T, typename ModelParameters, typename MomentumDistribution>
std::unique_ptr<MomentumDistribution> initialize_momentum_distribution(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, MomentumDistribution> &up_)
{
    return std::make_unique<MomentumDistribution>(std::complex<double>{up_.dist_sigma_real, up_.dist_sigma_imag}, up_.dist_sampler_n_autocorrelation, up_.dist_sampler_n_initialization, up_.dist_sampler_epsilon, 10100000);
}

template<typename T, typename ModelParameters>
std::unique_ptr<std::normal_distribution<double>> initialize_momentum_distribution(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, std::normal_distribution<double> > &up_)
{
    return std::make_unique<std::normal_distribution<double> >(0.0, 1.0);
}

template<typename T, typename ModelParameters>
std::unique_ptr<complex_gaussian_distribution_from_file> initialize_momentum_distribution(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, complex_gaussian_distribution_from_file > &up_)
{
    return std::make_unique<complex_gaussian_distribution_from_file>(up_.get_sigma(), 21000000, "ThrustGaussianDistribution");
}


template<typename T, typename ModelParameters, typename MomentumDistribution>
class SingleSiteCobridMonteCarloUpdate : public MCMCUpdateBase< SingleSiteCobridMonteCarloUpdate<T, ModelParameters, MomentumDistribution> >
{
public:
    explicit SingleSiteCobridMonteCarloUpdate(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, MomentumDistribution> &up_, typename ModelParameters::Model & model_) :
        up(up_), model(model_)
    {
        rand = std::uniform_real_distribution<double> (0,1);
        proposal_site = std::vector<T>{T(0)};
        momentum = std::vector<double>{0.0};
        momentum_backup = std::vector<double>{0.0};
        gaussian_dist_ptr = initialize_momentum_distribution(up_);
    }

    T operator() (const T site)
    {
        proposal_site[0] = site;
        auto current_energy = model.get_potential(site);

        // Sample momenta
        momentum[0] = (*gaussian_dist_ptr)(gen);
        momentum_backup[0] = momentum[0];

        // auto imag_energy = model.get_imag_potential(site)  + 0.5 * up.dist_sigma_imag * std::pow(momentum[0], 2.0) ;

        if(up.dist_sigma_imag != 1.0)
        {
            SingleSiteHamiltonianSystem hamiltonian_system(model);
            SingleSiteHamiltonianSystemMomentum hamiltonian_system_momentum(up.dist_sigma_imag);

            boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), std::make_pair(hamiltonian_system, hamiltonian_system_momentum),
                                                      make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                                                      0.0, up.dt, up.n);
        }
        else {
            SingleSiteHamiltonianSystem hamiltonian_system(model);

            boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                      make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                                                      0.0, up.dt, up.n);
        }

        // ToDo: It is possible to check a transition in between with the help of a observer

        proposal_site[0] = model.normalize(proposal_site[0]);
        auto proposal_energy = model.get_potential(proposal_site[0]);

        // auto imag_proposal_energy = model.get_imag_potential(proposal_site[0])  + 0.5 * up.dist_sigma_imag * std::pow(momentum[0], 2.0) ;
        // std::cout << imag_energy - imag_proposal_energy << std::endl;

        if(up.reweight_momentum == "true") {
            if (rand(gen) < std::min(
                    1.0, std::exp(-1.0 * (proposal_energy - current_energy) -
                        up.dist_sigma_real / 2.0 * (std::pow(momentum[0], 2.0) - std::pow(momentum_backup[0], 2.0)))))
                return proposal_site[0];
        }
        else
        {
            if (rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy - current_energy))))
                return proposal_site[0];
        }
        return site;
    }

    typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< std::vector<T> > symplectic_stepper;
private:

    struct SingleSiteHamiltonianSystem
    {
        SingleSiteHamiltonianSystem(typename ModelParameters::Model & model_) : hmc_model(model_)
        {}

        void operator()( const std::vector<T> &q , std::vector<T> &dpdt ) const {
            dpdt[0] = -1.0 * hmc_model.get_imag_drift_term(q[0]);
        }

        typename ModelParameters::Model & hmc_model;
    };

    struct SingleSiteHamiltonianSystemMomentum
    {
        SingleSiteHamiltonianSystemMomentum(const double imag_sigma_value_) : imag_sigma_value(imag_sigma_value_)
        {}

        void operator()( const std::vector<T> &p , std::vector<T> &dqdt ) const {
            dqdt[0] =  imag_sigma_value * p[0];
        }

        const double imag_sigma_value;
    };

    const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters, MomentumDistribution> & up;
    typename ModelParameters::Model & model;

    std::vector<double> proposal_site;
    std::vector<double> momentum;
    std::vector<double> momentum_backup;

    std::unique_ptr<MomentumDistribution> gaussian_dist_ptr;
    std::uniform_real_distribution<double> rand;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COBRID_MONTE_CARLO_HPP
