//
// Created by lukas on 07.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_HYBRID_MONTE_CARLO_HPP


#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../../mcmc_update_base.hpp"


template<typename T, typename ModelParameters>
class SingleSiteHybridMonteCarloUpdate;


template<typename T, typename ModelParameters>
class SingleSiteHybridMonteCarloUpdateParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteHybridMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                    dt(get_value_by_key<double>("dt", 0.01)),
                                                                    n(get_value_by_key<int>("n", 20))

    {}

    explicit SingleSiteHybridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : SingleSiteHybridMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteHybridMonteCarloUpdate";
    }

    typedef SingleSiteHybridMonteCarloUpdate<T, ModelParameters> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteHybridMonteCarloUpdate<T, ModelParameters>;

    const double dt;
    const int n;
};


template<typename T, typename ModelParameters>
class SingleSiteHybridMonteCarloUpdate : public MCMCUpdateBase< SingleSiteHybridMonteCarloUpdate<T, ModelParameters> >
{
public:
    explicit SingleSiteHybridMonteCarloUpdate(const SingleSiteHybridMonteCarloUpdateParameters<T, ModelParameters> &up_, typename ModelParameters::Model & model_) : up(up_), model(model_)
    {
        normal = std::normal_distribution<double> (0.0, 1.0);
        rand = std::uniform_real_distribution<double> (0,1);
        proposal_site = std::vector<T>{T(0)};
        momentum = std::vector<double>{0.0};
    }

    T operator() (const T site)
    {
        proposal_site[0] = site;

        // Sample momenta
        momentum[0] = normal(gen);

        SingleSiteHamiltonianSystem hamiltonian_system(model);

        boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                  make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                                                  0.0, up.dt, up.n);
        proposal_site[0] = model.normalize(proposal_site[0]);

        if(rand(gen) < std::min(1.0, std::exp(-1.0 * (model.get_potential(proposal_site[0]) - model.get_potential(site)))))
            return proposal_site[0];
        else
            return site;
    }

    T operator() (const T site, const std::vector< T* > neighbours)
    {
        proposal_site[0] = site;

        // Sample momenta
        momentum[0] = normal(gen);

        SiteHamiltonianSystemWithNeighbours hamiltonian_system(model, neighbours);

        boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                  make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                                                  0.0, up.dt, up.n);
        proposal_site[0] = model.normalize(proposal_site[0]);

        if(rand(gen) < std::min(1.0, std::exp(-1.0 * (model.get_potential(proposal_site[0], neighbours) - model.get_potential(site, neighbours)))))
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
            dpdt[0] = -1.0 * hmc_model.get_drift_term(q[0]);
        }

        typename ModelParameters::Model & hmc_model;
    };

    struct SiteHamiltonianSystemWithNeighbours
    {
        SiteHamiltonianSystemWithNeighbours(typename ModelParameters::Model & model_, const std::vector< T* > &neighbours_)
         : hmc_model(model_), neighbours(neighbours)
        {}

        void operator()( const std::vector<T> &q , std::vector<T> &dpdt ) const {
            dpdt[0] = -1.0 * hmc_model.get_drift_term(q[0], neighbours);
        }

        typename ModelParameters::Model & hmc_model;
        const std::vector< T* > &neighbours;
    };

    const SingleSiteHybridMonteCarloUpdateParameters<T, ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::vector<double> proposal_site;
    std::vector<double> momentum;
    std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> rand;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_HYBRID_MONTE_CARLO_HPP
