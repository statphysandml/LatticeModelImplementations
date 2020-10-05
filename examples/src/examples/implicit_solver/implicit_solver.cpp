#include "../../../include/examples/implicit_solver/implicit_solver.hpp"

template<typename ImplicitIntegralSolver>
void examples::implicit_solver::implicit_iterative_solve_simple_func()
{
    examples::implicit_solver::gaussian_function integrand;
    ImplicitIntegralSolver implict_integral_solver(integrand.lower_bound, integrand.upper_bound, 100, 10e-7);
    readdy::util::integration::integrator integrator;
    auto result = implict_integral_solver.solve(0.5, integrand, integrator);
    result = integrand.transformer(result);
    std::cout << "Upper integration boundary for an outcome of 0.5 for a gaussian integral: " << result << std::endl;
}


void examples::implicit_solver::run_implicit_integral_solvers()
{
    examples::implicit_solver::implicit_iterative_solve_simple_func<iterative_solver>();
    examples::implicit_solver::implicit_iterative_solve_simple_func<ceres_solver>();
}
