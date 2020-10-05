//
// Created by lukas on 17.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_THRUST_COMPLEX_GAUSSIAN_DISTRIBUTION_HPP
#define LATTICEMODELIMPLEMENTATIONS_THRUST_COMPLEX_GAUSSIAN_DISTRIBUTION_HPP

#ifdef THRUST

#include "thrust_header.hpp"

#include <thrust/random.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/complex.h>
#include <thrust/iterator/zip_iterator.h>
#include <curand_kernel.h>
#include <boost/filesystem.hpp>

#include <complex>
#include <random>

#include "param_helper/fileos.hpp"
#include "param_helper/pathfinder.hpp"

class thrust_complex_gaussian_distribution
{
public:
    thrust_complex_gaussian_distribution(std::complex<double> a_, uint n_autocorrelation_,
                                         uint n_initialization, const double epsilon_, const int M_, const std::string files_dir_="None");

    template< class Generator >
    double operator()( Generator& g )
    {
        return get_random_number();
    }

    std::complex<double> sample_complex()
    {
        return get_complex_random_number();
    }

    const static std::string name()
    {
        return "thrust_complex_gaussian_distribution";
    }

    unsigned long long get_total_updates() const
    {
        return total_updates;
    }

    void write_data_to_file(const std::string files_dir="ThrustComplexGaussian") const;

private:
    const std::complex<double> a;
    const uint n_autocorrelation;
    const int M;
    const double epsilon;
    const std::string files_dir;

    dev_vec real_random_numbers;
    dev_vec imag_random_numbers;
    // dev_vec_int indices;
    // thrust::permutation_iterator<dev_iterator, dev_iterator_int> iter_random_numbers;
    unsigned int i;

    thrust::host_vector<int> indices;
    thrust::permutation_iterator<thrust::host_vector<cudaT>::iterator, thrust::host_vector<int>::iterator> iter_real_random_numbers;
    thrust::permutation_iterator<thrust::host_vector<cudaT>::iterator, thrust::host_vector<int>::iterator> iter_imag_random_numbers;
    thrust::host_vector<cudaT> host_real_random_numbers;
    thrust::host_vector<cudaT> host_imag_random_numbers;

    unsigned long long total_updates;

    thrust::minstd_rand rng;
    thrust::uniform_real_distribution<float> dist;

    void update_random_numbers(const uint n_updates);
    double get_random_number();
    std::complex<double> get_complex_random_number();
    void evolve(const uint n_updates);
};

#endif

#endif //LATTICEMODELIMPLEMENTATIONS_THRUST_COMPLEX_GAUSSIAN_DISTRIBUTION_HPP
