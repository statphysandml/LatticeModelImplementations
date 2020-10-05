#include "../../../include/lattice_model_impl/thrust/thrust_finite_integration.hpp"

#ifdef THRUST

template<typename Func>
struct triangle
{
    Func f;
    const double lower_limit;
    const double interval;

    triangle(Func &f_, const double lower_limit_, const double interval_) : f(f_), lower_limit(lower_limit_), interval(interval_)
    {}

    __host__ __device__
    double operator() (const int& n)
    {
        return 0.5 * interval * (f(lower_limit + n * interval) + f(lower_limit + (n + 1) * interval));
    }
};


template<typename Func>
struct square
{
    Func f;
    const double lower_limit;
    const double interval;

    square(Func &f_, const double lower_limit_, const double interval_) : f(f_), lower_limit(lower_limit_), interval(interval_)
    {}

    __host__ __device__
    double operator() (const int& n)
    {
        return f(lower_limit + n * interval) * interval;
    }
};


inline readdy::util::thrust_finite_integration::integrator::integrator(const int n_) :
    n(n_)
{
    // results_sq.resize(n);
    results.resize(n);
}


template<typename Func, typename ScalarType>
inline std::pair<ScalarType, ScalarType> readdy::util::thrust_finite_integration::integrator::integrate(Func f, ScalarType lowerLimit, ScalarType upperLimit, bool watch)
{
    // std::cout << lowerLimit << " " << upperLimit << std::endl;

    const double interval = (upperLimit - lowerLimit) / n;

    // thrust::transform(thrust::make_counting_iterator(0), thrust::make_counting_iterator(n), results_sq.begin(), square<Func>(f, lowerLimit, interval));
    thrust::transform(thrust::make_counting_iterator(0), thrust::make_counting_iterator(n), results.begin(), triangle<Func>(f, lowerLimit, interval));

    // ScalarType integral_sq = thrust::reduce(results_sq.begin(), results_sq.end(), 0.0, thrust::plus<ScalarType>());
    ScalarType integral = thrust::reduce(results.begin(), results.end(), 0.0, thrust::plus<ScalarType>());

    // ScalarType integral_sq = thrust::transform_reduce(thrust::make_counting_iterator(0), thrust::make_counting_iterator(n), square<Func>(f, lowerLimit, interval), 0.0, thrust::plus<ScalarType>());
    // ScalarType integral = thrust::transform_reduce(thrust::make_counting_iterator(0), thrust::make_counting_iterator(n), triangle<Func>(f, lowerLimit, interval), 0.0, thrust::plus<ScalarType>());

    // std::cout << "Integrals" << integral << "\t" << integral_sq << std::endl;

    return std::make_pair(integral, 0.0); // abs(integral - integral_sq));
}
#endif