//
// Created by lukas on 30.09.20.
//

#ifndef COMPLEXMONTECARLO_COMPLEX_IMPLICIT_UPDATE_HPP
#define COMPLEXMONTECARLO_COMPLEX_IMPLICIT_UPDATE_HPP

#include <boost/math/tools/roots.hpp>

#include "../../mcmc_update_base.hpp"

template<typename Integrator, typename ImplicitIntegralSolver, typename TransitionRate, typename ModelParameters, typename SamplerCl>
class ComplexImplicitUpdate;

template<typename Integrator, typename ImplicitIntegralSolver, typename TransitionRate, typename ModelParameters, typename SamplerCl>
class ComplexImplicitUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit ComplexImplicitUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                   expansion(get_value_by_key<std::string>("expansion", "full")),
                                                                   n_grid_points(get_value_by_key<int>("n_grid_points", 10000))
    {}

    explicit ComplexImplicitUpdateParameters(const double eps_, const int n_grid_points_=10000)
        : ComplexImplicitUpdateParameters(json{{"n_grid_points", n_grid_points_}, {"eps", eps_}})
    {}

    static std::string name() {
        return "ComplexImplicitUpdate";
    }

    typedef ComplexImplicitUpdate<Integrator, ImplicitIntegralSolver, TransitionRate, ModelParameters, SamplerCl> MCMCUpdate;

private:
    friend class ComplexImplicitUpdate<Integrator, ImplicitIntegralSolver, TransitionRate, ModelParameters, SamplerCl>;

    const int n_grid_points;
    const std::string expansion;
};


template<typename Integrator, typename ImplicitIntegralSolver, typename TransitionRate, typename ModelParameters, typename SamplerCl>
class ComplexImplicitUpdate : public MCMCUpdateBase< ComplexImplicitUpdate<Integrator, ImplicitIntegralSolver, TransitionRate, ModelParameters, SamplerCl>, SamplerCl >
{
public:
    explicit ComplexImplicitUpdate(
            const ComplexImplicitUpdateParameters<Integrator, ImplicitIntegralSolver, TransitionRate, ModelParameters, SamplerCl> &up_,
            typename ModelParameters::Model & model_
    ) : MCMCUpdateBase<ComplexImplicitUpdate<Integrator, ImplicitIntegralSolver, TransitionRate, ModelParameters, SamplerCl>, SamplerCl>(up_.eps),
            up(up_), model(model_), integrator(Integrator(up.n_grid_points)), implicit_integral_solver(ImplicitIntegralSolver(-1.0, 1.0))
    {
        uniform = std::uniform_real_distribution<double>(0.0, 1.0);
    }

    template<typename T>
    T operator() (const T site)
    {
        // Compute based on site and based on proposal distribution the normalization factor - Check for monotonic increase
        // Find respective updated site
        TransitionRate integrand(model, this->sampler, site, 1.0, up.expansion);

        // Normalization
        auto normalization = compute_normalization_factor<TransitionRate>(integrand);
        integrand.set_normalization(normalization);

        implicit_integral_solver.update_integration_bounds(integrand.get_lower_bound(), integrand.get_upper_bound());

        // Implicit solver
        auto new_real_site = implicit_integral_solver.solve(uniform(gen), integrand, integrator);

        // Back transform
        new_real_site = integrand.transform(new_real_site);

        // Compute imaginary site
        double new_imag_site = integrand.get_new_imag_state({new_real_site, site.imag()});

        return std::complex<double>{new_real_site, new_imag_site};
    }

    double get_normalization_factor(const std::complex<double> site)
    {
        // Compute based on site and based on proposal distribution the normalization factor - Check for monotonic increase
        TransitionRate integrand(model, this->sampler, site, 1.0, up.expansion);
        return compute_normalization_factor<TransitionRate>(integrand);
    }

protected:
    const ComplexImplicitUpdateParameters<Integrator, ImplicitIntegralSolver, TransitionRate, ModelParameters, SamplerCl> & up;
    typename ModelParameters::Model & model;
    Integrator integrator;
    ImplicitIntegralSolver implicit_integral_solver;

    std::uniform_real_distribution<double> uniform;

    template<typename T>
    double compute_normalization_factor(TransitionRate &integrand)
    {
        auto normalization = integrand.compute_normalization_factor(integrator, integrand);
        if(normalization.second > 0.02)
        {
            auto state = integrand.get_state();
            std::cout << "Large error estimation in normalization factor detected: " << normalization.second
                      << " for state " << state.real() << " + i" << state.imag() << std::endl;
        }
        return normalization.first;
    }
};


#endif //COMPLEXMONTECARLO_COMPLEX_IMPLICIT_UPDATE_HPP
