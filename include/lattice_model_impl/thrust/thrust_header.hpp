//
// Created by lukas on 17.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_THRUST_HEADER_HPP
#define LATTICEMODELIMPLEMENTATIONS_THRUST_HEADER_HPP

#include <iostream>

#if defined(GPU) && defined(THRUST)

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace thrust::placeholders;

typedef double cudaT;

typedef thrust::device_vector<cudaT> dev_vec;
typedef thrust::device_vector<int> dev_vec_int;
typedef thrust::device_vector<bool> dev_vec_bool;
typedef thrust::host_vector< dev_vec_int * > dev_ptrvec_vec_int;

typedef thrust::device_vector<cudaT>::iterator dev_iterator;
typedef thrust::device_vector<int>::iterator dev_iterator_int;
typedef thrust::device_vector<bool>::iterator dev_iterator_bool;

typedef thrust::device_vector<cudaT>::const_iterator const_dev_iterator;
typedef thrust::device_vector<int>::const_iterator const_dev_iterator_int;
typedef thrust::device_vector<bool>::const_iterator const_dev_iterator_bool;

typedef thrust::device_vector< dev_iterator_int* > dev_iter_vec_int;

#elif THRUST

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace thrust::placeholders;

typedef double cudaT;

typedef thrust::host_vector<cudaT> dev_vec;
typedef thrust::host_vector<int> dev_vec_int;
typedef thrust::host_vector<bool> dev_vec_bool;
typedef thrust::host_vector< dev_vec_int * > dev_ptrvec_vec_int;

typedef thrust::host_vector<cudaT>::iterator dev_iterator;
typedef thrust::host_vector<int>::iterator dev_iterator_int;
typedef thrust::host_vector<bool>::iterator dev_iterator_bool;

typedef thrust::host_vector<cudaT>::const_iterator const_dev_iterator;
typedef thrust::host_vector<int>::const_iterator const_dev_iterator_int;
typedef thrust::host_vector<bool>::const_iterator const_dev_iterator_bool;

typedef thrust::host_vector< dev_iterator_int* > dev_iter_vec_int;

#endif


#ifdef THRUST

void print_system_info();

// from https://github.com/thrust/thrust/blob/master/examples/stream_compaction.cu

template <typename Iterator>
void print_range(const std::string& name, Iterator first, Iterator last)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;

    std::cout << name << ": ";
    thrust::copy(first, last, std::ostream_iterator<T>(std::cout, " "));
    std::cout << "\n";
}

template <typename Iterator>
void print_range_in_os(Iterator first, Iterator last, std::ofstream &os)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;
    thrust::copy(first, last, std::ostream_iterator<T>(os, " "));
}

#endif

#endif //LATTICEMODELIMPLEMENTATIONS_THRUST_HEADER_HPP
