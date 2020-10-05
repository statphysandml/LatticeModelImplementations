//
// Created by lukas on 21.02.20.
//

#ifndef MAIN_THRUST_FINITE_INTEGRATION_HPP
#define MAIN_THRUST_FINITE_INTEGRATION_HPP

#include <limits>
#include <numeric>

#ifdef THRUST

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/remove.h>
#include <thrust/complex.h>

#include "thrust_header.hpp"

namespace readdy {
    namespace util {
        namespace thrust_finite_integration {
            struct integrator
            {
                integrator(const int n_);

                template<typename Func, typename ScalarType>
                std::pair<ScalarType, ScalarType> integrate(Func f, ScalarType lowerLimit, ScalarType upperLimit, bool watch=false);

                const int n;
                dev_vec results;
            };
        }
    }
}

#endif

#endif //MAIN_THRUST_FINITE_INTEGRATION_HPP
