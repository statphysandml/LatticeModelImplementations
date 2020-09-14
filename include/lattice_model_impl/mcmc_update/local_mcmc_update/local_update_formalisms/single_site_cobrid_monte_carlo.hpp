//
// Created by lukas on 07.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COBRID_MONTE_CARLO_HPP


#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../../mcmc_update_base.hpp"

#include "../../../distribution/imaginary_gaussian_distribution.hpp"


template<typename T, typename ModelParameters>
class SingleSiteCobridMonteCarloUpdate;


template<typename T, typename ModelParameters>
class SingleSiteCobridMonteCarloUpdateParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteCobridMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                    dt(get_value_by_key<double>("dt", 0.01)),
                                                                    n(get_value_by_key<int>("n", 20)),
                                                                    dist_sigma_real(get_value_by_key<double>("dist_sigma_real", 0.001)),
                                                                    dist_sampler_n_autocorrelation(get_value_by_key<uint>("dist_sampler_n_autocorrelation", 1000)),
                                                                    dist_sampler_n_initialization(get_value_by_key<uint>("dist_sampler_n_initialization", 1000)),
                                                                    dist_sampler_epsilon(get_value_by_key<double>("dist_sampler_epsilon", 0.001))

    {}

    explicit SingleSiteCobridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : SingleSiteCobridMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteCobridMonteCarloUpdate";
    }

    typedef SingleSiteCobridMonteCarloUpdate<T, ModelParameters> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteCobridMonteCarloUpdate<T, ModelParameters>;

    const double dt;
    const int n;
    const double dist_sigma_real;
    const uint dist_sampler_n_autocorrelation;
    const uint dist_sampler_n_initialization;
    const double dist_sampler_epsilon;
};


template<typename T, typename ModelParameters>
class SingleSiteCobridMonteCarloUpdate : public MCMCUpdateBase< SingleSiteCobridMonteCarloUpdate<T, ModelParameters> >
{
public:
    explicit SingleSiteCobridMonteCarloUpdate(const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters> &up_, typename ModelParameters::Model & model_) :
        up(up_), model(model_), imag_gaussian_dist({up_.dist_sigma_real, 1.0}, {0.0, 0.0}, up_.dist_sampler_n_autocorrelation, up_.dist_sampler_n_initialization, up_.dist_sampler_epsilon)
    {
        rand = std::uniform_real_distribution<double> (0,1);
        proposal_site = std::vector<T>{T(0)};
        momentum = std::vector<double>{0.0};
        momentum_backup = std::vector<double>{0.0};
    }

    T operator() (const T site)
    {
        proposal_site[0] = site;
        auto current_energy = model.get_potential(site);

        // Sample momenta
        momentum[0] = imag_gaussian_dist(gen);
        momentum_backup[0] = momentum[0];

        SingleSiteHamiltonianSystem hamiltonian_system(model);

        boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                  make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                                                  0.0, up.dt, up.n);

        proposal_site[0] = model.normalize(proposal_site[0]);
        auto proposal_energy = model.get_potential(proposal_site[0]);

        if(rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy - current_energy))))
            return proposal_site[0];
        else
            return site;
    }

    T operator() (const T site, const std::vector< T* > neighbours)
    {
        proposal_site[0] = site;
        auto current_energy = model.get_potential(site, neighbours);

        // Sample momenta
        momentum[0] = imag_gaussian_dist(gen);
        momentum_backup[0] = momentum[0];

        SiteHamiltonianSystemWithNeighbours hamiltonian_system(model, neighbours);

        boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                  make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                                                  0.0, up.dt, up.n);
        proposal_site[0] = model.normalize(proposal_site[0]);
        auto proposal_energy = model.get_potential(proposal_site[0], neighbours);

        if(rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy - current_energy))))
            return proposal_site[0];
        else
            return site;
    }

    typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< std::vector<T> > symplectic_stepper;
protected:

    struct SingleSiteHamiltonianSystem
    {
        SingleSiteHamiltonianSystem(typename ModelParameters::Model & model_) : hmc_model(model_)
        {}

        void operator()( const std::vector<T> &q , std::vector<T> &dpdt ) const {
            dpdt[0] = -1.0 * hmc_model.get_imag_drift_term(q[0]);
        }

        typename ModelParameters::Model & hmc_model;
    };

    struct SiteHamiltonianSystemWithNeighbours
    {
        SiteHamiltonianSystemWithNeighbours(typename ModelParameters::Model & model_, const std::vector< T* > &neighbours_)
                : hmc_model(model_), neighbours(neighbours)
        {}

        void operator()( const std::vector<T> &q , std::vector<T> &dpdt ) const {
            dpdt[0] = -1.0 * hmc_model.get_imag_drift_term(q[0], neighbours);
        }

        typename ModelParameters::Model & hmc_model;
        const std::vector< T* > &neighbours;
    };

    const SingleSiteCobridMonteCarloUpdateParameters<T, ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::vector<double> proposal_site;
    std::vector<double> momentum;
    std::vector<double> momentum_backup;
    imaginary_gaussian_distribution imag_gaussian_dist;
    std::uniform_real_distribution<double> rand;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COBRID_MONTE_CARLO_HPP
