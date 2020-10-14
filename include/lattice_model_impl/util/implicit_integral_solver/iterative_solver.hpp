//
// Created by lukas on 29.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_ITERATIVE_SOLVER_HPP
#define LATTICEMODELIMPLEMENTATIONS_ITERATIVE_SOLVER_HPP

#include <iostream>

#include "mcmc_simulation/util/random.hpp"

struct iterative_solver
{
    iterative_solver(const double lower_bound_, const double upper_bound_, const double maximum_number_iterations_ = 100, const double accuracy_=10e-5)
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
        double current_upper_bound = upper_bound;
        double current_lower_bound = lower_bound;

        double current_guess = (current_upper_bound + current_lower_bound) / 2.0;
        auto outcome = integrator.integrate(func, lower_bound, current_guess, true);

        double integral = outcome.first;
        // double uncertainty = outcome.second;

        auto c = 0;
        while (true) {
            if (current_upper_bound - current_lower_bound < accuracy)
                break;
            if (c > 30) {
                std::cout << "c > 30 " << c << std::endl;
            }
            if (c > maximum_number_iterations) {
                std::cout << "Implicit equation did not converge!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (integral < target) {
                current_lower_bound = std::max(lower_bound, current_guess - 0.3 * accuracy); // - uncertainty);
            } else {
                current_upper_bound = std::min(upper_bound, current_guess + 0.3 * accuracy); // + uncertainty);
            }
            current_guess = (current_upper_bound + current_lower_bound) / 2.0;

            outcome = integrator.integrate(func, lower_bound, current_guess, true);
            integral = outcome.first;
            // uncertainty = outcome.second;
            c++;
        }

        if(abs(outcome.first - target) > 0.01)
            std::cout << "Large difference between target value and integral: " << outcome.first - target << " with target " << target << " and current guess " << current_guess << std::endl;

        return current_guess;
    }

    double lower_bound;
    double upper_bound;

    const double maximum_number_iterations;
    const double accuracy;
};

#endif //LATTICEMODELIMPLEMENTATIONS_ITERATIVE_SOLVER_HPP
