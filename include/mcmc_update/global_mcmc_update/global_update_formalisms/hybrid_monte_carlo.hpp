//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP

#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../global_update.hpp"


template<typename T, typename ModelParameters>
class HybridMonteCarloUpdate;


template<typename T, typename ModelParameters>
class HybridMonteCarloUpdateParameters : public GlobalUpdateFormalismParameters {
public:
    explicit HybridMonteCarloUpdateParameters(const json params_) : GlobalUpdateFormalismParameters(params_),
                                                                    dt(get_value_by_key<double>("dt", 0.01)),
                                                                    n(get_value_by_key<int>("n", 20))

    {}

    explicit HybridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : HybridMonteCarloUpdateParameters(json{
        {"dt", dt_},
        {"n", n_}})
    {}

    static std::string name() {
        return "HybridMonteCarloUpdate";
    }

    typedef HybridMonteCarloUpdate<T, ModelParameters> UpdateFormalism;

protected:
    friend class HybridMonteCarloUpdate<T, ModelParameters>;

    const double dt;
    const int n;
};


template<typename T, typename ModelParameters>
class HybridMonteCarloUpdate : public GlobalUpdateFormalism< HybridMonteCarloUpdate<T, ModelParameters> >
{
public:
    explicit HybridMonteCarloUpdate(const HybridMonteCarloUpdateParameters<T, ModelParameters> &up_, typename ModelParameters::Model & model_) : up(up_), model(model_)
    {
        normal = std::normal_distribution<double> (0.0, 1.0);
        rand = std::uniform_real_distribution<double> (0,1);
    }

    template<typename Lattice>
    void initialize_global_update(const Lattice& lattice)
    {
        momenta = std::vector<double> (lattice.size(), 0.0);
    }

    template<typename Lattice>
    void operator() (Lattice& lattice)
    {
        auto current_energy = lattice.energy();
        auto current_lattice_grid(lattice.get_system_representation());

        typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< std::vector<T> > symplectic_stepper;

        // Sample momenta
        std::generate(momenta.begin(), momenta.end(), [this] () { return normal(gen); });

        HamiltonianSystem hamiltonian_system( model, lattice.get_neighbours());

        boost::numeric::odeint::integrate_n_steps( symplectic_stepper() , hamiltonian_system ,
                         make_pair( boost::ref( lattice.get_system_representation() ) , boost::ref( momenta ) ) ,
                         0.0 , up.dt, up.n);

        lattice.normalize(lattice.get_system_representation());

        auto proposal_energy = lattice.energy();
        if(rand(gen) >= std::min(1.0, std::exp(-1.0 * (proposal_energy - current_energy))))
        {
            // ToDo: Rewrite?
            auto& lattice_grid = lattice.get_system_representation();
            lattice_grid = current_lattice_grid;
        }
    }

protected:

    struct HamiltonianSystem
    {
        HamiltonianSystem(typename ModelParameters::Model & model_, const std::vector< std::vector < T* > > &neighbours_) : hmc_model(model_), hmc_neighbours(neighbours_)
        {}

        void operator()( const std::vector<T> &q , std::vector<T> &dpdt ) const {
            for(auto i=0; i < q.size(); i++)
            {
                dpdt[i] = -1.0 * hmc_model.get_drift_term(q[i], hmc_neighbours[i]);
            }
        }

        typename ModelParameters::Model & hmc_model;
        const std::vector< std::vector < T* > > &hmc_neighbours;
    };

    struct SingleSiteHamiltonianSystem
    {
        SingleSiteHamiltonianSystem(typename ModelParameters::Model & model_) : hmc_model(model_)
        {}

        void operator()( const std::vector<T> &q , std::vector<T> &dpdt ) const {
            dpdt[0] = -1.0 * hmc_model.get_drift_term(q[0]);
        }

        typename ModelParameters::Model & hmc_model;
    };

    const HybridMonteCarloUpdateParameters<T, ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::vector<double> momenta;
    std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> rand;
};

#endif //LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
