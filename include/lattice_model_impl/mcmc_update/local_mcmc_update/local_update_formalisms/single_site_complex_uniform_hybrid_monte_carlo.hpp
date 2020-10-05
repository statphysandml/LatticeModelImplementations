//
// Created by lukas on 02.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_UNIFORM_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_UNIFORM_HYBRID_MONTE_CARLO_HPP


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

#include "single_site_complex_hybrid_monte_carlo.hpp"

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexUniformHybridMonteCarlo;

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexUniformHybridMonteCarloParameters;

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexUniformHybridMonteCarloParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteComplexUniformHybridMonteCarloParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                                     dt(get_value_by_key<double>("dt", 0.01)),
                                                                                     n(get_value_by_key<int>("n", 20)),
                                                                                     dist_sigma_real(get_value_by_key<double>("dist_sigma_real", 0.001)),
                                                                                     dist_sigma_imag(get_value_by_key<double>("dist_sigma_imag", 1.0)),
                                                                                     dist_sampler_n_autocorrelation(get_value_by_key<uint>("dist_sampler_n_autocorrelation", 10000)),
                                                                                     dist_sampler_n_initialization(get_value_by_key<uint>("dist_sampler_n_initialization", 10000)),
                                                                                     dist_sampler_epsilon(get_value_by_key<double>("dist_sampler_epsilon", 0.0001))
    {}

    explicit SingleSiteComplexUniformHybridMonteCarloParameters(const double dt_, const int n_
    ) : SingleSiteComplexUniformHybridMonteCarloParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteComplexUniformHybridMonteCarlo";
    }

    std::complex<double> get_sigma() const
    {
        return {dist_sigma_real, dist_sigma_imag};
    }

    typedef SingleSiteComplexUniformHybridMonteCarlo<ModelParameters, MomentumDistribution> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteComplexUniformHybridMonteCarlo<ModelParameters, MomentumDistribution>;
    friend struct momentum_distribution<ModelParameters, MomentumDistribution>;

    const double dt;
    const int n;
    const double dist_sigma_real;
    const double dist_sigma_imag;
    const uint dist_sampler_n_autocorrelation;
    const uint dist_sampler_n_initialization;
    const double dist_sampler_epsilon;
};

template<typename ModelParameters, typename MomentumDistribution>
class SingleSiteComplexUniformHybridMonteCarlo : public MCMCUpdateBase< SingleSiteComplexUniformHybridMonteCarlo<ModelParameters, MomentumDistribution> >
{
public:
    explicit SingleSiteComplexUniformHybridMonteCarlo(const SingleSiteComplexUniformHybridMonteCarloParameters<ModelParameters, MomentumDistribution> &up_, typename ModelParameters::Model & model_) :
            up(up_), model(model_), gaussian_dist(up)
    {
        /* proposal_site = std::vector<double>{0.0, 0.0, 0.0, 0.0};
        current_site = std::vector<double>{0.0, 0.0, 0.0, 0.0}; */
        proposal_site = std::vector<double>{0.0, 0.0};
        current_site = std::vector<double>{0.0, 0.0};
    }

    std::complex<double> operator() (const std::complex<double> site)
    {
        // std::cout << "Site " << site.real() << " + i" << site.imag() << std::endl;

        // X
        proposal_site[0] = site.real();
        // Sample momenta
        auto sample = gaussian_dist.sample_complex(gen);
        proposal_site[1] = sample.real();

        current_site = proposal_site;

        // Preparations for real update

        hamilton_system_H_imag_constant hamiltonian_system(model, site.imag(), sample.imag(), up.dist_sigma_real, up.dist_sigma_imag);

        auto energy = hamiltonian_system.get_energies(proposal_site);
        double imag_energy = hamiltonian_system.get_H_imag_energy(proposal_site);

        boost::numeric::odeint::integrate_n_steps(stepper_type(),
                                                  hamiltonian_system,
                                                  proposal_site,
                                                  0.0, up.dt, up.n);

        auto proposal_energy = hamiltonian_system.get_energies(proposal_site);
        double proposal_imag_energy = hamiltonian_system.get_H_imag_energy(proposal_site);

        // Preparations for imag update
        std::vector<double> imag_proposal_site(2, 0);
        imag_proposal_site[0] = site.real();
        sample = gaussian_dist.sample_complex(gen);
        imag_proposal_site[1] = sample.real();

        hamilton_system_H_imag_constant hamiltonian_system_imag(model, site.imag(), sample.imag(), up.dist_sigma_real, up.dist_sigma_imag);

        boost::numeric::odeint::integrate_n_steps(stepper_type(),
                                                  hamiltonian_system_imag,
                                                  imag_proposal_site,
                                                  0.0, up.dt, up.n);

        auto new_imag_state = hamiltonian_system.compute_proposal_imag_state(proposal_site[0], imag_proposal_site[0], proposal_energy.second, energy.second);

        // std::cout << proposal_energy.first - energy.first << " + i" << proposal_imag_energy - imag_energy << std::endl;

        if(std::abs(proposal_imag_energy -imag_energy) > 0.1 or std::abs(proposal_energy.first - energy.first) > 120)
            std::cout << "Something is going wrong?" << std::endl;

        return hamiltonian_system.accept_reject_step(proposal_site[0], current_site[0], new_imag_state, proposal_energy.first, energy.first);
    }

    typedef boost::numeric::odeint::runge_kutta4< std::vector<double> > stepper_type;
private:

    struct hamilton_system_H_imag_constant
    {
        hamilton_system_H_imag_constant(typename ModelParameters::Model & model_, const double current_imag_,
                                        const double current_imag_momentum_, const double real_sigma_, const double imag_sigma_) :
                                        hmc_model(model_), real_sigma(real_sigma_), imag_sigma(imag_sigma_),
                                        current_imag(current_imag_), current_imag_momentum(current_imag_momentum_)
        {
            rand = std::uniform_real_distribution<double> (0,1);
        }

        void operator()( const std::vector<double> &x , std::vector<double> &dxdt, const double /* t */  ) const {
            auto drift = hmc_model.get_drift_term(std::complex<double> {x[0], current_imag});
            dxdt[0] = real_sigma * current_imag_momentum + imag_sigma * x[1];
            dxdt[1] = -1.0 * drift.imag();
        }

        std::pair<double, double> get_energies(std::vector<double> &x) const
        {
            auto energy = hmc_model.get_potential(std::complex<double>{x[0], current_imag});
            return std::pair<double, double> (energy.real() + (0.5 * std::complex<double> (real_sigma, imag_sigma) * std::pow(std::complex<double> {x[1], current_imag_momentum}, 2.0)).real(),
                    energy.imag());
        }

        double get_H_imag_energy(std::vector<double> &x) const
        {
            auto energy = hmc_model.get_potential(std::complex<double>{x[0], current_imag});
            return energy.imag() + (0.5 * std::complex<double> (real_sigma, imag_sigma) * std::pow(std::complex<double> {x[1], current_imag_momentum}, 2.0)).imag();
        }

        std::complex<double> accept_reject_step(double new_site, double site, double new_imag_state, double proposal_energy, double energy)
        {
            if(rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy - energy)))) //  and std::abs(proposal_energy.imag() - current_energy.imag()) < up.max_delta_s_imag)
                return {new_site, new_imag_state};
            else
                return {site, new_imag_state};
        }

        double compute_proposal_imag_state(double proposal_state, double imag_proposal_site, double imag_proposal_energy, double imag_energy)
        {
            return current_imag - (imag_proposal_site - proposal_state) * tan(-0.5 * (imag_proposal_energy - imag_energy));
        }

        typename ModelParameters::Model & hmc_model;
        const double current_imag;
        const double current_imag_momentum;
        const double real_sigma;
        const double imag_sigma;
        std::uniform_real_distribution<double> rand;
    };

    const SingleSiteComplexUniformHybridMonteCarloParameters<ModelParameters, MomentumDistribution> & up;
    typename ModelParameters::Model & model;

    std::vector<double> proposal_site;
    std::vector<double> current_site;

    momentum_distribution<ModelParameters, MomentumDistribution> gaussian_dist;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_UNIFORM_HYBRID_MONTE_CARLO_HPP
