#include "../../../include/examples/integration/integration.hpp"

#include <complex>

struct test_func
{
    double operator() (const double x)
    {
        return std::exp(std::complex<double>{(std::pow(x, 2) * std::sin(x) * std::exp(-1.0) * std::sin(x) * std::exp(-1.0) * std::sin(x) * std::exp(-1.0)), 0.0}).real();
    }
};

template<typename Integrator>
void test_integration(const int n=100)
{
    test_func thrust_func;

    Integrator integrator(n);

    std::pair<double, double> result;
    for(auto i = 0; i < 100000; i++) {
        result = integrator.integrate(thrust_func, -1.3, 3.3);
    }
    std::cout << result.first << " " << result.second << std::endl;
}


void integrate()
{
    std::cout << "Integration test" << std::endl;

    test_integration<readdy::util::integration::integrator>();
}