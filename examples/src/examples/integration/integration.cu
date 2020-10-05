#include "../../../include/examples/integration/integration.hpp"

#include "../../../../src/lattice_model_impl/thrust/thrust_integration.cu"
#include "../../../../src/lattice_model_impl/thrust/thrust_finite_integration.cu"

struct test_func : public thrust::unary_function<double, double>
{
    __host__ __device__
    double operator() (const double x)
    {
        return thrust::exp(thrust::complex<double>{(std::pow(x, 2) * std::sin(x) * std::exp(-1.0) * std::sin(x) * std::exp(-1.0) * std::sin(x) * std::exp(-1.0)), 0.0}).real();
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
    test_integration<readdy::util::thrust_integration::integrator<double>>();
    test_integration<readdy::util::thrust_finite_integration::integrator>();
}