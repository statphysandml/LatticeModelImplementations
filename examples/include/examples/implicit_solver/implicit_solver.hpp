//
// Created by lukas on 30.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONSEXAMPLES_IMPLICIT_SOLVER_HPP
#define LATTICEMODELIMPLEMENTATIONSEXAMPLES_IMPLICIT_SOLVER_HPP

#include "lattice_model_impl/thrust/thrust_header.hpp"
#include "lattice_model_impl/util/implicit_integral_solver_header.hpp"
#include "lattice_model_impl/util/integration_header.hpp"

namespace examples {
    namespace implicit_solver {

        struct gaussian_function
        {
            #ifdef THRUST
            __host__ __device__
            #endif
            double operator() (const double x)
            {
                return jacobian(x) * 1/std::sqrt(2 * M_PI) * std::exp(-0.5 * std::pow(transformer(x), 2.0));
            }

            #ifdef THRUST
            __host__ __device__
            #endif
            double jacobian(const double x)
            {
                return -2.0 / (std::pow(x, 2.0) - 1.0);
            }

            struct transformer_func
            {
                #ifdef THRUST
                __host__ __device__
                #endif
                double operator() (const double val)
                {
                    return std::log((1.0 + val) / (1.0 - val));
                }
            };

            transformer_func transformer;

            const double lower_bound = -1.0;
            const double upper_bound = 1.0;
        };

        template<typename ImplicitIntegralSolver>
        void implicit_iterative_solve_simple_func();
        void run_implicit_integral_solvers();

        #ifdef THRUST
        template<typename ImplicitIntegralSolver>
        void thrust_implicit_iterative_solve_simple_func();
        void thrust_run_implicit_integral_solvers();
        #endif
    }
}

#endif //LATTICEMODELIMPLEMENTATIONSEXAMPLES_IMPLICIT_SOLVER_HPP
