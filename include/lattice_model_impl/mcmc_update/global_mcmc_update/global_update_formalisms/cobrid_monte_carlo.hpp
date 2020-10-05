//
// Created by lukas on 02.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_COBRID_MONTE_CARLO_HPP

#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../../mcmc_update_base.hpp"

#include "../../../util/distribution/imaginary_gaussian_distribution.hpp"


template<typename T, typename ModelParameters>
class CobridMonteCarloUpdate;


template<typename T, typename ModelParameters>
class CobridMonteCarloUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit CobridMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                        dt(get_value_by_key<double>("dt")),
                                                                        n(get_value_by_key<int>("n")),
                                                                        dist_sigma_real(get_value_by_key<double>("dist_sigma_real", 0.001)),
                                                                        dist_sampler_n_autocorrelation(get_value_by_key<uint>("dist_sampler_n_autocorrelation", 1000)),
                                                                        dist_sampler_n_initialization(get_value_by_key<uint>("dist_sampler_n_initialization", 1000)),
                                                                        dist_sampler_epsilon(get_value_by_key<double>("dist_sampler_epsilon", 0.001))


    {}

    explicit CobridMonteCarloUpdateParameters(const double dt_, const int n_
    ) : CobridMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "CobridMonteCarloUpdate";
    }

    typedef CobridMonteCarloUpdate<T, ModelParameters> MCMCUpdate;

protected:
    friend class CobridMonteCarloUpdate<T, ModelParameters>;

    const double dt;
    const int n;
    const double dist_sigma_real;
    const uint dist_sampler_n_autocorrelation;
    const uint dist_sampler_n_initialization;
    const double dist_sampler_epsilon;
};


template<typename T, typename ModelParameters>
class CobridMonteCarloUpdate : public MCMCUpdateBase< CobridMonteCarloUpdate<T, ModelParameters> >
{
public:
    explicit CobridMonteCarloUpdate(const CobridMonteCarloUpdateParameters<T, ModelParameters> &up_, typename ModelParameters::Model & model_) :
            up(up_), model(model_), imag_gaussian_dist({up_.dist_sigma_real, 1.0}, up_.dist_sampler_n_autocorrelation, up_.dist_sampler_n_initialization, up_.dist_sampler_epsilon)
    {
        rand = std::uniform_real_distribution<double> (0,1);
    }

    template<typename Lattice>
    void initialize_mcmc_update(const Lattice& lattice)
    {
        momenta = std::vector<double> (lattice.size(), 0.0);
    }

    template<typename Lattice>
    void operator() (Lattice& lattice)
    {
        auto current_lattice_grid(lattice.get_system_representation());
        auto current_energy = lattice.energy();

        typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< std::vector<T> > symplectic_stepper;

        // Sample momenta
        std::generate(momenta.begin(), momenta.end(), [this] () { return imag_gaussian_dist(gen); }); // imaginary_gaussian_distribution(gen)
        std::copy(momenta.begin(), momenta.end(), std::back_inserter(momenta_backup));

        HamiltonianSystem hamiltonian_system(model, lattice.get_neighbours());

        // std::cout << std::pow(momenta[0], 2.0) / 2.0 + model.get_imag_potential(lattice.get_system_representation());
        boost::numeric::odeint::integrate_n_steps(symplectic_stepper(), hamiltonian_system,
                                                  make_pair(boost::ref( lattice.get_system_representation() ), boost::ref(momenta)),
                                                  0.0, up.dt, up.n);

        lattice.normalize(lattice.get_system_representation());
        // std::cout << " == " << std::pow(momenta[0], 2.0) / 2.0 + model.get_imag_potential(q[0]) << std::endl;

        auto proposal_energy = lattice.energy();
        auto& site = lattice.get_system_representation();
        if(rand(gen) >= std::min(1.0,
                                 std::exp(-1.0 * (proposal_energy - current_energy)))) // - 1.0 / 2.0 * (
            // std::pow(momenta[0], 2.0) - std::pow(momenta_backup[0], 2.0)))))
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
                dpdt[i] = -1.0 * hmc_model.get_imag_drift_term(q[i], hmc_neighbours[i]);
            }
        }

        typename ModelParameters::Model & hmc_model;
        const std::vector< std::vector < T* > > &hmc_neighbours;
    };

    const CobridMonteCarloUpdateParameters<T, ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::vector<double> momenta;
    std::vector<double> momenta_backup;
    imaginary_gaussian_distribution imag_gaussian_dist;
    std::uniform_real_distribution<double> rand;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COBRID_MONTE_CARLO_HPP
