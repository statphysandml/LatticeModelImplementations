#include "../../../include/lattice_model_impl/thrust/thrust_integration.hpp"

#ifdef THRUST

template<typename Func>
struct transfo
{
    Func f;
    const double midpoint;
    const double halfInterval;

    transfo(Func f_, const double midpoint_, const double halfInterval_) : f(f_), midpoint(midpoint_), halfInterval(halfInterval_)
    {}

    __host__ __device__
    double operator() (const double& abscissa, const double& weight)
    {
        return halfInterval * f(halfInterval * abscissa + midpoint) * weight;
    }
};

template<typename ScalarType>
struct negate
{
    __host__ __device__
    ScalarType operator() (const ScalarType& val)
    {
        return -1.0 * val;
    }
};

template<typename scalar_type>
inline readdy::util::thrust_integration::integrator<scalar_type>::integrator(const int) :
    len(abscissaeGaussKronrod201<scalar_type>.size()),lenGauss(abscissaeGauss201<scalar_type>.size())
{
    thrustweightsGaussKronrod201.resize(2 * len - 1);
    thrustabscissaeGaussKronrod201.resize(2 * len - 1);
    thrustweightsGauss201.resize(2 * lenGauss);

    resultsGauss.resize(2 * lenGauss);
    results.resize(2 * len - 1);
    abscissae.resize(2 * len - 1);
    indices.resize(2 * len - 1);
    thrust::sequence(indices.begin(), indices.end());

    thrust::copy(weightsGauss201<scalar_type>.begin(), weightsGauss201<scalar_type>.end(), thrustweightsGauss201.begin());
    thrust::copy(weightsGauss201<scalar_type>.rbegin(), weightsGauss201<scalar_type>.rend(), thrustweightsGauss201.begin() + lenGauss);

    thrust::copy(weightsGaussKronrod201<scalar_type>.begin(), weightsGaussKronrod201<scalar_type>.end(), thrustweightsGaussKronrod201.begin());
    thrust::copy(weightsGaussKronrod201<scalar_type>.rbegin() + 1, weightsGaussKronrod201<scalar_type>.rend(), thrustweightsGaussKronrod201.begin() + len);

    thrust::copy(abscissaeGaussKronrod201<scalar_type>.begin(), abscissaeGaussKronrod201<scalar_type>.end(), thrustabscissaeGaussKronrod201.begin());

    // Set weights
    thrust::transform(thrustabscissaeGaussKronrod201.rbegin() + len, thrustabscissaeGaussKronrod201.rend(), thrustabscissaeGaussKronrod201.begin() + len, negate<scalar_type>());
}

template<typename scalar_type>
template<typename Func, typename ScalarType>
inline std::pair<ScalarType, ScalarType> readdy::util::thrust_integration::integrator<scalar_type>::integrate(Func f, ScalarType lowerLimit, ScalarType upperLimit, bool watch)
{
    if (lowerLimit > upperLimit) {
        throw std::invalid_argument("lower limit cannot be larger than upper limit");
    } else if (lowerLimit == upperLimit) {
        return std::make_pair(static_cast<ScalarType>(0.), static_cast<ScalarType>(0.));
    }

    const ScalarType midpoint = (lowerLimit + upperLimit) / 2.;
    const ScalarType halfInterval = (upperLimit - lowerLimit) / 2.;

/*    thrust::transform(iw.thrustabscissaeGaussKronrod201.begin(), iw.thrustabscissaeGaussKronrod201.end(), iw.abscissae.begin(),
            [halfInterval, midpoint] __host__ __device__ (const ScalarType &abscissa)
            {
                return halfInterval * abscissa + midpoint;
            });

    // thrust::transform(iw.abscissae.begin(), iw.abscissae.end(), iw.abscissae.begin(), f);*/
    thrust::transform(thrustabscissaeGaussKronrod201.begin(), thrustabscissaeGaussKronrod201.end(), thrustweightsGaussKronrod201.begin(), results.begin(), transfo<Func>(f, midpoint, halfInterval));
    /* thrust::transform(iw.abscissae.begin(), iw.abscissae.end(), iw.thrustweightsGaussKronrod201.begin(), iw.results.begin(),
                      [halfInterval] __host__ __device__ (const ScalarType& val, const ScalarType& weight) { return halfInterval * val * weight; }); */
    ScalarType integralKronrod = thrust::reduce(results.begin(), results.end(), 0.0, thrust::plus<ScalarType>());

    /* auto values_end = thrust::remove_copy_if(iw.abscissae.begin(), iw.abscissae.end(), iw.indices.begin(), iw.abscissae.begin(), is_even<int>());
    iw.abscissae.resize(values_end - iw.abscissae.begin());

    thrust::transform(iw.abscissae.begin(), iw.abscissae.end(), iw.thrustweightsGauss201.begin(), iw.resultsGauss.begin(),
                      [halfInterval] __host__ __device__ (const ScalarType& val, const ScalarType& weight) { return halfInterval * val * weight; });
    ScalarType integralGauss = thrust::reduce(iw.resultsGauss.begin(), iw.resultsGauss.end(), 0.0, thrust::plus<ScalarType>()); */

    const ScalarType absoluteErrorEstimate = std::abs(integralKronrod - 1.0);
    return std::make_pair(integralKronrod, 0.0); // absoluteErrorEstimate);
}

// template class readdy::util::thrust_integration::IntegrationWeights<double>;

#endif