//
// Created by lukas on 29.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_CERES_SOLVER_HPP
#define LATTICEMODELIMPLEMENTATIONS_CERES_SOLVER_HPP

#include <iostream>
#include "ceres/ceres.h"

#include "mcmc_simulation/util/random.hpp"

template<typename T, typename Func, typename Integrator>
struct CeresNumericCostFunctor {
    explicit CeresNumericCostFunctor(const T& target_, Func &func_, Integrator &integrator_, const double& lower_bound_) :
            target(target_), integrator(integrator_), func(func_), lower_bound(lower_bound_)
    {}

    bool operator()(const T* const x, T* residual) const {
        auto outcome = integrator.integrate(func, lower_bound, x[0], true);
        double integral = outcome.first;
        // double uncertainty = outcome.second;

        residual[0] = integral - target;
        return true;
    }

    const T& target;
    Func& func;
    Integrator& integrator;
    const double &lower_bound;
};


struct ceres_solver
{
    ceres_solver(const double lower_bound_, const double upper_bound_, const double maximum_number_iterations_ = 100, const double accuracy_=10e-5)
            : lower_bound(lower_bound_), upper_bound(upper_bound_), maximum_number_iterations(maximum_number_iterations_), accuracy(accuracy_)
    {}

    void update_integration_bounds(double lower_bound_, double upper_bound_)
    {
        upper_bound = upper_bound_;
        lower_bound = lower_bound_;
    }

    template<typename T, typename Func, typename Integrator>
    T solve(const T& target, Func &func, Integrator &integrator)
    {
        double current_guess = (upper_bound + lower_bound) / 2.0;

        ceres::Problem problem;
        ceres::CostFunction* cost_function =
                new ceres::NumericDiffCostFunction<CeresNumericCostFunctor< T, Func, Integrator >, ceres::CENTRAL, 1, 1>(
                new CeresNumericCostFunctor< T, Func, Integrator >(target, func, integrator, lower_bound));

        problem.AddResidualBlock(cost_function, NULL, &current_guess);

        // Run the solver!
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        // options.minimizer_progress_to_stdout = true;
        ceres::Solver::Summary summary;
        Solve(options, &problem, &summary);

        auto outcome = integrator.integrate(func, lower_bound, current_guess, true);

        if(abs(outcome.first - target) > 0.001)
            std::cout << "Large difference between target value and integral: " << outcome.first - target << " with target " << target << " and current guess " << current_guess << std::endl;

        return current_guess;
    }

    double lower_bound;
    double upper_bound;

    const double maximum_number_iterations;
    const double accuracy;
};

#endif //LATTICEMODELIMPLEMENTATIONS_CERES_SOLVER_HPP
