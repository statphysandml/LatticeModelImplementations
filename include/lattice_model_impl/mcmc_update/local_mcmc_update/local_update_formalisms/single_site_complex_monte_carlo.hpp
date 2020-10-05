//
// Created by lukas on 18.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_MONTE_CARLO_HPP

#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
// #include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"
// #include "boost/numeric/odeint/integrate/integrate_adaptive.hpp"
#include "boost/numeric/odeint/stepper/generation.hpp"
// #include "boost/numeric/odeint.hpp"
// #include "boost/numeric/odeint/"

#include "../../mcmc_update_base.hpp"


template<typename ModelParameters, typename SamplerCl>
class SingleSiteComplexMonteCarloUpdate;


template<typename ModelParameters, typename SamplerCl>
class SingleSiteComplexMonteCarloUpdateParameters : public MCMCUpdateBaseParameters  {
public:
    explicit SingleSiteComplexMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                              dt(get_value_by_key<double>("dt", 0.01)),
                                                                              n(get_value_by_key<int>("n", 20)),
                                                                              epsilon(get_value_by_key<double>("epsilon", 0.001)),
                                                                              max_delta_s_imag(get_value_by_key<double>("max_delta_s_imag", 0.01))

    {}

    explicit SingleSiteComplexMonteCarloUpdateParameters(const double dt_, const int n_
    ) : SingleSiteComplexMonteCarloUpdateParameters(json{
            {"dt", dt_},
            {"n", n_}})
    {}

    static std::string name() {
        return "SingleSiteComplexMonteCarloUpdate";
    }

    typedef SingleSiteComplexMonteCarloUpdate<ModelParameters, SamplerCl> MCMCUpdate;
    typedef typename ModelParameters::Model Model;

protected:
    friend class SingleSiteComplexMonteCarloUpdate<ModelParameters, SamplerCl>;

    const double dt;
    const int n;
    const double epsilon;
    const double max_delta_s_imag;
};


template<typename ModelParameters, typename SamplerCl>
class SingleSiteComplexMonteCarloUpdate : public MCMCUpdateBase< SingleSiteComplexMonteCarloUpdate<ModelParameters, SamplerCl> >
{
public:
    explicit SingleSiteComplexMonteCarloUpdate(const SingleSiteComplexMonteCarloUpdateParameters<ModelParameters, SamplerCl> &up_, typename ModelParameters::Model & model_) :
            up(up_), model(model_)
    {
        rand = std::uniform_real_distribution<double> (0,1);
        normal = std::normal_distribution<double> (0, 1);
        rand_sign = std::uniform_int_distribution<int>(0, 1);
        rand_n = std::uniform_int_distribution<int>(1, up.n);
        // q = std::vector<double>{0.0};
        // p = std::vector<double>{0.0};
    }

    std::complex<double> operator() (const std::complex<double> site)
    {
        auto current_energy = model.get_potential(site);

        // q[0] = site.real();
        // p[0] = site.imag();

        auto proposal_site = std::complex<double> {site.real(), site.imag()};

        proposal_site.real(proposal_site.real() + up.epsilon * normal(gen));

        double sign = 1.0;
        // double sign = 1.0;
        // auto n = rand_n(gen);

        auto imag_energy = current_energy.imag();

        /* boost::numeric::odeint::integrate_adaptive(
                boost::numeric::odeint::make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-6 ) ,
                                                   hamilton_system(model, sign) , proposal_site , 0.0 , 10.0 , 0.01 ); */

        boost::numeric::odeint::integrate_n_steps(stepper_type(),
                                                  hamilton_system(model, sign),
                                                  proposal_site,
                                                  0.0, up.dt, up.n);

        proposal_site.real(proposal_site.real() + up.epsilon * normal(gen));

        auto proposal_energy = model.get_potential(proposal_site);

        auto imag_proposal_energy = proposal_energy.imag();
        std::cout << imag_energy - imag_proposal_energy << std::endl;

        if(rand(gen) < std::min(1.0, std::exp(-1.0 * (proposal_energy.real() - current_energy.real()))) and std::abs(proposal_energy.imag() - current_energy.imag()) < up.max_delta_s_imag)
            return proposal_site;
        else
            return site;
    }
    // typedef boost::numeric::odeint::runge_kutta_cash_karp54< std::complex<double> > error_stepper_type;
    typedef boost::numeric::odeint::runge_kutta4< std::complex<double> > stepper_type;
protected:

    struct hamilton_system
    {
        hamilton_system(typename ModelParameters::Model & model_, const double sign_) : hmc_model(model_), sign(sign_)
        {}

        void operator()( const std::complex<double> &x , std::complex<double> &dxdt, const double /* t */  ) const {
            dxdt.real(sign * hmc_model.get_imag_const_imag_drift_term(x));
            dxdt.imag(-1.0 * sign * hmc_model.get_real_const_imag_drift_term(x));
        }

        const double sign;
        typename ModelParameters::Model & hmc_model;
    };

    /* struct hamilton_system_x
    {
        explicit hamilton_system_x(typename ModelParameters::Model & model_, const double sign_) : hmc_model(model_), sign(sign_)
        {}

        void operator()( const std::vector<double> &p_ , std::vector<double> &dqdt ) const {
            dqdt[0] = sign * hmc_model.get_imag_const_imag_drift_term(p_[0]);
        }

        const double sign;
        typename ModelParameters::Model & hmc_model;
    }; */

    /* struct hamilton_system_y
    {
        explicit hamilton_system_y(typename ModelParameters::Model & model_, const double sign_) : hmc_model(model_), sign(sign_)
        {}

        void operator()( const std::vector<double> &q_ , std::vector<double> &dpdt ) const {
            dpdt[0] = -1.0 * sign * hmc_model.get_real_const_imag_drift_term(q_[0]);
        }

        const double sign;
        typename ModelParameters::Model & hmc_model;
    }; */

    const SingleSiteComplexMonteCarloUpdateParameters<ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;

    std::uniform_real_distribution<double> rand;
    std::normal_distribution<double> normal;
    std::uniform_int_distribution<int> rand_sign;
    std::uniform_int_distribution<int> rand_n;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SINGLE_SITE_COMPLEX_MONTE_CARLO_HPP
