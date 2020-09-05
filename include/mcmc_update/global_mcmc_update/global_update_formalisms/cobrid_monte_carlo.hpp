//
// Created by lukas on 02.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_COBRID_MONTE_CARLO_HPP

#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../global_update.hpp"

#include "../../../distribution/imaginary_gaussian_distribution.hpp"


template<typename T, typename ModelParameters>
class CobridMonteCarloUpdate;


template<typename T, typename ModelParameters>
class CobridMonteCarloUpdateParameters : public GlobalUpdateFormalismParameters {
public:
    explicit CobridMonteCarloUpdateParameters(const json params_) : GlobalUpdateFormalismParameters(params_),
                                                                    dt(get_value_by_key<double>("dt")),
                                                                    n(get_value_by_key<int>("n"))

    {}

    explicit CobridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : CobridMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "CobridMonteCarloUpdate";
    }

    typedef CobridMonteCarloUpdate<T, ModelParameters> UpdateFormalism;

protected:
    friend class CobridMonteCarloUpdate<T, ModelParameters>;

    const double dt;
    const int n;
};


template<typename T, typename ModelParameters>
class CobridMonteCarloUpdate : public GlobalUpdateFormalism< CobridMonteCarloUpdate<T, ModelParameters> >
{
public:
    explicit CobridMonteCarloUpdate(const CobridMonteCarloUpdateParameters<T, ModelParameters> &up_, typename ModelParameters::Model & model_) :
        up(up_), model(model_)
    {
        // normal = std::normal_distribution<double> (0.0, 1.0);
        rand = std::uniform_real_distribution<double> (0,1);
    }

    template<typename Lattice>
    void initialize_global_update(const Lattice& lattice)
    {
        momenta = std::vector<double> (lattice.size(), 0.0);

        if(Lattice::get_name() == "Lattice") // ToDo: Change to enum
            hamilton_system_type = "coupled";
        else
            hamilton_system_type = "single_site";
    }

    template<typename Lattice>
    void operator() (Lattice& lattice)
    {
        auto current_lattice_grid(lattice.get_system_representation());
        auto current_energy = model.get_potential(current_lattice_grid);

        typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< std::vector<T> > symplectic_stepper;

        // Sample momenta
        std::generate(momenta.begin(), momenta.end(), [this] () { return normal(gen); });
        /* if(hamilton_system_type == "coupled") {
            CoupledHamiltonianSystem hamiltonian_system(model, lattice.get_neighbours());

            boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                      make_pair(boost::ref( lattice.get_system_representation() ), boost::ref(momenta)),
                                                      0.0, up.dt, up.n);
        }
        else
        { */
        // std::cout << std::pow(momenta[0], 2.0) / 2.0 + model.get_potential(lattice.get_system_representation());

            SingleSiteHamiltonianSystem hamiltonian_system(model);

            std::vector<T> q {lattice.get_system_representation()};
            boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                      make_pair(boost::ref(q), boost::ref(momenta)),
                                                      0.0, up.dt, up.n);
            lattice.normalize(q[0]);
        // }
        // std::cout << " == " << std::pow(momenta[0], 2.0) / 2.0 + model.get_potential(q[0]) << std::endl;
        auto& site = lattice.get_system_representation();
        site = q[0];
        auto proposal_energy = model.get_potential(site);
        if(rand(gen) >= std::min(1.0, std::exp(-1.0 * (proposal_energy - current_energy))))
        {
            // ToDo: Rewrite?
            site = current_lattice_grid;
        }
    }

protected:

    struct CoupledHamiltonianSystem
    {
        CoupledHamiltonianSystem(typename ModelParameters::Model & model_, const std::vector< std::vector < T* > > &neighbours_) : hmc_model(model_), hmc_neighbours(neighbours_)
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

    const CobridMonteCarloUpdateParameters<T, ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::vector<double> momenta;
    std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> rand;

    std::string hamilton_system_type;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COBRID_MONTE_CARLO_HPP
