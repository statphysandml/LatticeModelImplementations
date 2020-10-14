//
// Created by lukas on 08.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_HYBRID_MONTE_CARLO_PURE_REAL_SAMPLING_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_HYBRID_MONTE_CARLO_PURE_REAL_SAMPLING_HPP


#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../../mcmc_update_base.hpp"


// HMC without momentum, just keep the imaginary part of the action constant


template<typename ModelParameters, typename SamplerCl>
class SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate;


template<typename ModelParameters, typename SamplerCl>
class SingleSiteComplexHybridMonteCarloPureRealSamplingUpdateParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteComplexHybridMonteCarloPureRealSamplingUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                               dt(get_value_by_key<double>("dt", 0.01)),
                                                                               n(get_value_by_key<int>("n", 20))

    {}

    explicit SingleSiteComplexHybridMonteCarloPureRealSamplingUpdateParameters(const double dt_, const int n_
    ) : SingleSiteComplexHybridMonteCarloPureRealSamplingUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate";
    }

    typedef SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate<ModelParameters, SamplerCl> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate<ModelParameters, SamplerCl>;

    const double dt;
    const int n;
};


template<typename ModelParameters, typename SamplerCl>
class SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate : public MCMCUpdateBase< SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate<ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate(const SingleSiteComplexHybridMonteCarloPureRealSamplingUpdateParameters<ModelParameters, SamplerCl> &up_, typename ModelParameters::Model & model_)
        : MCMCUpdateBase< SingleSiteComplexHybridMonteCarloPureRealSamplingUpdate<ModelParameters, SamplerCl>, SamplerCl >(up_.eps), up(up_), model(model_)
    {
        rand = std::uniform_real_distribution<double> (0,1);
        normal = std::normal_distribution<double> (0, 1);
        rand_n = std::uniform_int_distribution<int>(1, up.n);

        proposal_site = std::vector<double>{0.0};
        proposal_site_imag_update = std::vector<double>{0.0};
        momentum = std::vector<double>{0.0};
    }

    std::complex<double> operator() (const std::complex<double> site)
    {
        auto current_energy = model.get_potential(site);

        // Real part
        proposal_site[0] = site.real();
        momentum[0] = normal(gen);

        auto hmc_current_energy = current_energy.imag() + 0.5 * std::pow(momentum[0], 2.0);

        boost::numeric::odeint::integrate_n_steps(
                symplectic_stepper(), hamilton_system(model, site.imag()),
                make_pair(boost::ref(proposal_site), boost::ref(momentum)),
                0.0, up.dt, up.n);

        auto proposal_energy = model.get_potential({proposal_site[0], site.imag()});

        auto hmc_proposal_energy = proposal_energy.imag() + 0.5 * std::pow(momentum[0], 2.0);

        std::cout << "Energy difference: " << hmc_current_energy - hmc_proposal_energy << std::endl;

        // Imaginary part
        proposal_site_imag_update[0] = site.real();
        momentum[0] = normal(gen);

        boost::numeric::odeint::integrate_n_steps(
                symplectic_stepper(), hamilton_system(model, site.imag()),
                make_pair(boost::ref(proposal_site_imag_update), boost::ref(momentum)),
                0.0, up.dt, up.n);

        auto new_imag_site = compute_new_imag_state(site.imag(), proposal_energy.imag(), current_energy.imag());

        if(rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy.real() - current_energy.real()))))
            return {proposal_site[0], new_imag_site};
        else
            return {site.real(), new_imag_site};
    }

    typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< std::vector<double> > symplectic_stepper;
private:

    struct hamilton_system
    {
        hamilton_system(typename ModelParameters::Model & model_, const double current_imag_site_)
            : hmc_model(model_), current_imag_state(current_imag_site_)
        {}

        void operator()( const std::vector<double> &x , std::vector<double> &dxdt) const {
            dxdt[0] = -1.0 * hmc_model.get_drift_term(std::complex<double>{x[0], current_imag_state}).imag();
        }

        typename ModelParameters::Model & hmc_model;
        const double current_imag_state;
    };

    double compute_new_imag_state(const double current_imag, double imag_proposal_energy, double imag_energy)
    {
        // Does ist really have to be "-proposal_site[0]" or is "-current[0]" correct?
        return current_imag - (proposal_site_imag_update[0] - proposal_site[0]) * tan(-0.5 * (imag_proposal_energy - imag_energy));
    }

    const SingleSiteComplexHybridMonteCarloPureRealSamplingUpdateParameters<ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;

    std::vector<double> proposal_site;
    std::vector<double> proposal_site_imag_update;
    std::vector<double> momentum;

    std::uniform_real_distribution<double> rand;
    std::normal_distribution<double> normal;
    std::uniform_int_distribution<int> rand_n;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_HYBRID_MONTE_CARLO_PURE_REAL_SAMPLING_HPP
